"""
CURE (Clinical-failure Understanding & Redesign Engine) — Agent 0~5 통합 파이프라인
Colab에서 순서대로 검증한 로직을 스크립트 하나로 정리한 버전.

필요 패키지:
    pip install chembl_webresource_client admet-ai mmpdb rdkit requests pandas numpy
"""

import json
import os
import subprocess
import time

import numpy as np
import pandas as pd
import requests
from admet_ai import ADMETModel
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Contrib.SA_Score import sascorer

# ── 공용 클라이언트 / 모델 준비 ──────────────────────────────
molecule = new_client.molecule
mechanism = new_client.mechanism
target = new_client.target
activity = new_client.activity
similarity = new_client.similarity

model = ADMETModel()
remover = SaltRemover()

params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog(params)

TOXICITY_ENDPOINTS = [
    "AMES", "DILI", "hERG", "ClinTox", "Carcinogens_Lagunin", "Skin_Reaction",
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
    "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
]

# Deep-PK 항목명 매핑 — 확인된 것만 등록 (필요시 추가)
DEEPPK_KEY_MAP = {
    "DILI": "[Toxicity/Liver Injury I (DILI)] Predictions",
    "hERG": "[Toxicity/hERG Blockers] Predictions",
}

# mmpdb index의 --max-variable-heavies 기본값(10)을 그대로 쓰던 걸 명시적으로 드러냄.
# 이 값보다 원자 수가 많은 "변형 부위(variable fragment)"를 가진 매칭쌍은 index 단계에서
# 조용히 제외된다 — 실측 확인: Cisapride 사례에서 10(기본값)일 땐 후보 17개, 제한을 풀면
# (--max-variable-heavies none) 130개로 7.6배 늘어남. 즉 지금 이 값은 "국소적인 치환"만
# 후보로 보겠다는 실질적인 설계 선택이다. 더 큰 구조 변경까지 후보로 보고 싶다면 이 값을
# 올리거나 None으로 설정(제한 없음)할 것 — 다만 후보가 늘어난 만큼 mmpdb 빌드/transform
# 시간도 늘어난다.
MMPDB_MAX_VARIABLE_HEAVIES = 10


# ── 공용 유틸 ────────────────────────────────────────────────
RETRY_DELAYS = (2, 4, 8, 16)  # 초 단위 지수 백오프 (총 5회 시도)


def _with_retry(fn, label, delays=RETRY_DELAYS):
    # ChEMBL/ADMET-AI/Deep-PK 외부 서버 500/404/타임아웃 대응.
    # label에 화합물/단계 정보를 담아 실패 시 어디서 죽었는지 로그로 남긴다.
    last_exc = None
    for attempt, delay in enumerate((0,) + delays, start=1):
        if delay:
            print(f"  [재시도 {attempt - 1}/{len(delays)}] {label} — {delay}초 대기 후 재시도")
            time.sleep(delay)
        try:
            return fn()
        except Exception as e:
            last_exc = e
            print(f"  [에러] {label} 실패 (시도 {attempt}/{len(delays) + 1}): {type(e).__name__}: {e}")
    raise RuntimeError(f"{label} — {len(delays) + 1}회 시도 모두 실패: {last_exc}") from last_exc


def _work_path(compound_key, filename):
    # mmpdb 산출물은 화합물별 디렉터리로 분리해 배치 실행 시 파일 충돌을 방지한다.
    d = os.path.join("mmpdb_work", str(compound_key))
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, filename)


def _run_mmpdb(args, label):
    try:
        subprocess.run(args, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"  [mmpdb 실패] {label}\n    cmd: {' '.join(args)}\n    stderr: {e.stderr}")
        raise


def desalt(smiles):
    # 대이온(counter-ion) 제거
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit이 파싱할 수 없는 SMILES (ChEMBL 응답 확인 필요): {smiles!r}")
    return Chem.MolToSmiles(remover.StripMol(mol))


def strip_isotopes(smiles):
    # 동위원소 표지 제거 (구조 정규화용)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit이 파싱할 수 없는 SMILES (ChEMBL 응답 확인 필요): {smiles!r}")
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolToSmiles(mol)


def check_structural_alerts(smiles):
    # PAINS/Brenk 구조 경보 목록 반환
    mol = Chem.MolFromSmiles(smiles)
    return [m.GetDescription() for m in catalog.GetMatches(mol)]


def get_sa_score(smiles):
    # 합성 가능성 점수 (1=쉬움 ~ 10=어려움)
    mol = Chem.MolFromSmiles(smiles)
    return sascorer.calculateScore(mol)


def poll_deeppk(job_id, delay=5):
    # 진행중(dict)/완료(이중 인코딩 문자열) 둘 다 처리
    while True:
        r = _with_retry(
            lambda: requests.get(
                "https://biosig.lab.uq.edu.au/deeppk/api/predict",
                files={"job_id": (None, job_id)},
                timeout=30,
            ),
            label=f"Deep-PK 폴링(job_id={job_id})",
        )
        parsed = r.json()
        if isinstance(parsed, dict) and parsed.get("status") == "running":
            time.sleep(delay)
            continue
        return json.loads(parsed) if isinstance(parsed, str) else parsed


def deeppk_predict(smiles, key):
    resp = _with_retry(
        lambda: requests.post(
            "https://biosig.lab.uq.edu.au/deeppk/api/predict",
            files={"smiles": (None, smiles), "pred_type": (None, "toxicity")},
            timeout=30,
        ),
        label=f"Deep-PK 예측 요청(smiles={smiles[:30]}...)",
    )
    result = poll_deeppk(resp.json()["job_id"])
    if not isinstance(result, dict) or "0" not in result or key not in result["0"]:
        actual = list(result.get("0", {}).keys()) if isinstance(result, dict) and "0" in result else result
        raise KeyError(f"Deep-PK 응답에 예상 키({key!r}) 없음. 실제 응답 구조: {actual}")
    return result["0"][key]


# ── Agent 0: 타겟 자동탐색 ───────────────────────────────────
def get_chembl_id(drug_name):
    hits = molecule.search(drug_name)
    first = _with_retry(lambda: hits[0] if hits else None, label=f"ChEMBL ID 조회: {drug_name}")
    return first["molecule_chembl_id"] if first else None


def get_target_from_mechanism(chembl_id):
    hits = mechanism.filter(molecule_chembl_id=chembl_id)
    return _with_retry(lambda: hits[0] if hits else None, label=f"타겟 메커니즘 조회: {chembl_id}")


def get_target_name(target_chembl_id):
    return _with_retry(
        lambda: target.get(target_chembl_id)["pref_name"], label=f"타겟명 조회: {target_chembl_id}"
    )


def get_smiles(chembl_id):
    raw = _with_retry(
        lambda: molecule.get(chembl_id)["molecule_structures"]["canonical_smiles"],
        label=f"SMILES 조회: {chembl_id}",
    )
    return desalt(raw)


# ── Agent 1: 진단 (물성 / 효능 / 독성) ───────────────────────
def get_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "logP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "MW": Descriptors.MolWt(mol),
    }


def check_lipinski(props):
    violated = props["logP"] > 5 or props["MW"] > 500 or props["TPSA"] >= 140
    return "물성_문제" if violated else "정상"


def diagnose_properties(smiles):
    props = get_properties(smiles)
    if props is None:
        return {"물성_판정": "SMILES_파싱_실패"}
    props["물성_판정"] = check_lipinski(props)
    return props


def _parse_activity_values(resolved, context):
    # standard_value가 하나만 비정상(비수치 문자열 등)이어도 전체가 죽지 않도록 그 레코드만 건너뛴다.
    values = []
    for h in resolved:
        if h["standard_value"] is None:
            continue
        try:
            values.append(float(h["standard_value"]))
        except (TypeError, ValueError):
            print(f"  [건너뜀] 활성값 파싱 실패({context}): standard_value={h['standard_value']!r}")
    return values


def get_own_activity(molecule_chembl_id, target_chembl_id):
    hits = activity.filter(
        molecule_chembl_id=molecule_chembl_id,
        target_chembl_id=target_chembl_id,
        standard_type="IC50",
        standard_units="nM",
        pchembl_value__isnull=False,
    ).only(["standard_value"])
    resolved = _with_retry(
        lambda: list(hits), label=f"자체 활성값 조회: {molecule_chembl_id}/{target_chembl_id}"
    )
    values = _parse_activity_values(resolved, f"{molecule_chembl_id}/{target_chembl_id}")
    return min(values) if values else None


def get_target_activities(target_chembl_id, limit=2000):
    hits = activity.filter(
        target_chembl_id=target_chembl_id,
        standard_type="IC50",
        standard_units="nM",
        pchembl_value__isnull=False,
    ).only(["standard_value"])[:limit]
    resolved = _with_retry(lambda: list(hits), label=f"타겟 활성값 분포 조회: {target_chembl_id}")
    return _parse_activity_values(resolved, target_chembl_id)


def potency_percentile(own_value, population_values):
    if own_value is None or len(population_values) < 5:
        return None
    return float((np.array(population_values) < own_value).mean() * 100)


def diagnose_efficacy(molecule_chembl_id, target_chembl_id):
    own_value = get_own_activity(molecule_chembl_id, target_chembl_id)
    population = get_target_activities(target_chembl_id)
    pct = potency_percentile(own_value, population)
    if pct is None:
        return {"효능_판정": "데이터_부족"}
    verdict = "효력_부족" if pct >= 75 else "정상"
    return {"백분위": round(pct, 1), "효능_판정": verdict}


def diagnose_toxicity(preds):
    def risk_level(score):
        if score > 0.7:
            return "고위험"
        elif score >= 0.3:
            return "중위험"
        return "저위험"

    risks = {k: risk_level(preds[k]) for k in TOXICITY_ENDPOINTS}
    verdict = "독성_문제" if "고위험" in risks.values() else "정상"
    return {"항목별_위험도": risks, "독성_판정": verdict}


def get_primary_toxicity_endpoint(preds):
    # 가장 심각한(점수 높은) 독성 항목 하나 선택, 0.7 이하면 문제 없음
    scores = {k: preds[k] for k in TOXICITY_ENDPOINTS}
    worst = max(scores, key=scores.get)
    return worst if scores[worst] > 0.7 else None


def run_agent0_1(drug_name):
    # 타겟 미확정 시 효능 축만 스킵, 물성·독성 축은 계속 진행
    chembl_id = get_chembl_id(drug_name)
    if chembl_id is None:
        return None
    mech = get_target_from_mechanism(chembl_id)
    target_id = mech["target_chembl_id"] if mech else None
    smiles = get_smiles(chembl_id)
    preds = _with_retry(lambda: model.predict(smiles=smiles), label=f"ADMET-AI 예측: {drug_name}")
    return {
        "chembl_id": chembl_id,
        "target_id": target_id,
        "smiles": smiles,
        "물성": diagnose_properties(smiles),
        "효능": diagnose_efficacy(chembl_id, target_id) if target_id else {"효능_판정": "타겟_판단불가"},
        "독성": diagnose_toxicity(preds),
        "property_name": get_primary_toxicity_endpoint(preds),
    }


# ── Agent 2: 문제 부위 탐지 (mmpdb 빌드) ─────────────────────
def get_similar_compounds(smiles, threshold=70):
    hits = similarity.filter(smiles=smiles, similarity=threshold).only(
        ["molecule_chembl_id", "similarity"]
    )
    return _with_retry(lambda: list(hits), label=f"유사 화합물 검색(threshold={threshold})")


def get_similar_smiles(similar_compounds):
    # 유사 화합물 중 일부가 RDKit이 못 읽는 SMILES를 갖고 있어도(ChEMBL 데이터 품질 이슈)
    # 그 하나 때문에 화합물 전체가 죽지 않도록 건너뛰고 계속 진행한다.
    result = []
    for c in similar_compounds:
        cid = c["molecule_chembl_id"]
        raw = _with_retry(
            lambda cid=cid: molecule.get(cid)["molecule_structures"]["canonical_smiles"],
            label=f"유사 화합물 SMILES 조회: {cid}",
        )
        try:
            result.append((desalt(raw), cid))
        except ValueError as e:
            print(f"  [건너뜀] 유사 화합물 {cid}: {e}")
    return result


def dedup_against_query(pairs, query_smiles):
    seen = {strip_isotopes(query_smiles)}
    result = []
    for smi, cid in pairs:
        norm = strip_isotopes(smi)
        if norm not in seen:
            seen.add(norm)
            result.append((smi, cid))
    return result


def run_agent2(smiles, chembl_id, property_name):
    # mmpdb 산출물은 화합물(chembl_id)별 디렉터리에 써서 배치 실행 시 파일 충돌을 피한다.
    compounds_path = _work_path(chembl_id, "compounds.smi")
    fragdb_path = _work_path(chembl_id, "fragments.fragdb")
    mmpdb_path = _work_path(chembl_id, "mmp.db")
    props_path = _work_path(chembl_id, "props.txt")

    similar = get_similar_compounds(smiles, threshold=70)
    pairs = get_similar_smiles(similar)
    deduped = dedup_against_query(pairs, smiles)
    all_compounds = [(smiles, chembl_id)] + deduped

    with open(compounds_path, "w") as f:
        for smi, cid in all_compounds:
            f.write(f"{smi}\t{cid}\n")
    _run_mmpdb(
        ["mmpdb", "fragment", compounds_path, "-o", fragdb_path],
        label=f"fragment({chembl_id})",
    )
    max_var_heavies = "none" if MMPDB_MAX_VARIABLE_HEAVIES is None else str(MMPDB_MAX_VARIABLE_HEAVIES)
    _run_mmpdb(
        ["mmpdb", "index", fragdb_path, "--max-variable-heavies", max_var_heavies, "-o", mmpdb_path],
        label=f"index({chembl_id})",
    )

    batch_preds = _with_retry(
        lambda: model.predict(smiles=[s for s, _ in all_compounds]),
        label=f"ADMET-AI 배치 예측: {chembl_id} 비교군 {len(all_compounds)}개",
    )
    with open(props_path, "w") as f:
        f.write(f"ID\t{property_name}\n")
        for (s, cid), val in zip(all_compounds, batch_preds[property_name]):
            f.write(f"{cid}\t{val}\n")
    _run_mmpdb(
        ["mmpdb", "loadprops", "-p", props_path, mmpdb_path],
        label=f"loadprops({chembl_id})",
    )


# ── Agent 3: 구조 변형 제안 ───────────────────────────────────
def run_agent3(smiles, chembl_id, property_name):
    mmpdb_path = _work_path(chembl_id, "mmp.db")
    results_path = _work_path(chembl_id, "transform_results.csv")
    _run_mmpdb(
        ["mmpdb", "transform", "--smiles", smiles, mmpdb_path,
         "--property", property_name, "-o", results_path],
        label=f"transform({chembl_id})",
    )
    df = pd.read_csv(results_path, sep="\t")
    avg_col = f"{property_name}_avg"
    count_col = f"{property_name}_count"
    return df.sort_values(avg_col)[["SMILES", avg_col, count_col]].reset_index(drop=True)


# ── Agent 4 + 5: 재평가 게이트 + 랭킹 ─────────────────────────
def run_agent4_5(smiles, variants, property_name):
    avg_col = f"{property_name}_avg"
    baseline_alerts = set(check_structural_alerts(smiles))
    variants = variants.copy()
    variants["신규_구조_경보"] = variants["SMILES"].apply(
        lambda s: [a for a in check_structural_alerts(s) if a not in baseline_alerts]
    )
    variants["SAScore"] = variants["SMILES"].apply(get_sa_score)

    if property_name in DEEPPK_KEY_MAP:
        key = DEEPPK_KEY_MAP[property_name]
        variants["deeppk_판정"] = variants["SMILES"].apply(lambda s: deeppk_predict(s, key))
        variants["Agent4_통과"] = variants.apply(
            lambda r: len(r["신규_구조_경보"]) == 0 and r["deeppk_판정"] == "Safe",
            axis=1,
        )
    else:
        variants["Agent4_통과"] = variants["신규_구조_경보"].apply(len).eq(0)

    return variants[variants["Agent4_통과"]].sort_values(avg_col).reset_index(drop=True)


# ── 마스터 파이프라인 ─────────────────────────────────────────
def run_full_pipeline(drug_name, top_n_steps=(5, 15, 30)):
    try:
        diag = run_agent0_1(drug_name)
    except Exception as e:
        raise RuntimeError(f"[{drug_name}] Agent0-1(진단) 단계 실패: {e}") from e
    if diag is None:
        return {"status": "타겟_판단불가"}
    if diag["property_name"] is None:
        return {"status": "독성_문제_없음", "진단": diag}

    try:
        run_agent2(diag["smiles"], diag["chembl_id"], diag["property_name"])
    except Exception as e:
        raise RuntimeError(f"[{drug_name}] Agent2(MMP DB 구축) 단계 실패: {e}") from e

    try:
        all_variants = run_agent3(diag["smiles"], diag["chembl_id"], diag["property_name"])
    except Exception as e:
        raise RuntimeError(f"[{drug_name}] Agent3(구조 변형 제안) 단계 실패: {e}") from e

    checked = 0
    for n in top_n_steps:
        new_batch = all_variants.iloc[checked:n]
        checked = n
        if len(new_batch) == 0:
            continue
        try:
            ranked = run_agent4_5(diag["smiles"], new_batch, diag["property_name"])
        except Exception as e:
            raise RuntimeError(
                f"[{drug_name}] Agent4-5(재검증+랭킹) 단계 실패 (후보 {checked}개 확인 중): {e}"
            ) from e
        if len(ranked) > 0:
            return {
                "status": "완료",
                "타겟문제": diag["property_name"],
                "랭킹": ranked,
                "확인한_후보수": n,
            }

    return {"status": "통과_후보_없음", "타겟문제": diag["property_name"], "확인한_후보수": checked}


def run_batch(drug_names):
    # 여러 화합물 순차 실행(파일 충돌 방지를 위해 병렬화하지 않음).
    # 하나 에러 나도 중단 없이 기록만 남기고 계속 진행.
    results = {}
    for name in drug_names:
        print(f"=== {name} 실행 시작 ===")
        try:
            results[name] = run_full_pipeline(name)
        except Exception as e:
            print(f"[배치 에러] {name}: {e}")
            results[name] = {"status": "에러", "message": str(e)}
    return results


if __name__ == "__main__":
    # 검증됐던 3개(Cisapride, Troglitazone, Bromfenac) + hERG/DILI 사례 + 안전 대조군
    batch_results = run_batch([
        "Cisapride", "Troglitazone", "Bromfenac",
        "Terfenadine", "Astemizole",       # hERG 계열
        "Trovafloxacin",                   # DILI 계열
        "Metformin", "Amoxicillin", "Loratadine", "Ibuprofen",  # 안전 대조군
    ])
    for name, r in batch_results.items():
        print(name, "->", r.get("status"), "| 문제축:", r.get("타겟문제"))
