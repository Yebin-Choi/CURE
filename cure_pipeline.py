"""
CURE (Clinical-failure Understanding & Redesign Engine) — Agent 0~5 통합 파이프라인
Colab에서 순서대로 검증한 로직을 스크립트 하나로 정리한 버전.

필요 패키지:
    pip install chembl_webresource_client admet-ai mmpdb rdkit requests pandas numpy
"""

import json
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


# ── 공용 유틸 ────────────────────────────────────────────────
def desalt(smiles):
    # 대이온(counter-ion) 제거
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(remover.StripMol(mol))


def strip_isotopes(smiles):
    # 동위원소 표지 제거 (구조 정규화용)
    mol = Chem.MolFromSmiles(smiles)
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
        r = requests.get(
            "https://biosig.lab.uq.edu.au/deeppk/api/predict",
            files={"job_id": (None, job_id)},
        )
        parsed = r.json()
        if isinstance(parsed, dict) and parsed.get("status") == "running":
            time.sleep(delay)
            continue
        return json.loads(parsed) if isinstance(parsed, str) else parsed


def deeppk_predict(smiles, key):
    resp = requests.post(
        "https://biosig.lab.uq.edu.au/deeppk/api/predict",
        files={"smiles": (None, smiles), "pred_type": (None, "toxicity")},
    )
    return poll_deeppk(resp.json()["job_id"])["0"][key]


# ── Agent 0: 타겟 자동탐색 ───────────────────────────────────
def get_chembl_id(drug_name):
    hits = molecule.search(drug_name)
    return hits[0]["molecule_chembl_id"] if hits else None


def get_target_from_mechanism(chembl_id):
    hits = mechanism.filter(molecule_chembl_id=chembl_id)
    return hits[0] if hits else None


def get_target_name(target_chembl_id):
    return target.get(target_chembl_id)["pref_name"]


def get_smiles(chembl_id):
    raw = molecule.get(chembl_id)["molecule_structures"]["canonical_smiles"]
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


def get_own_activity(molecule_chembl_id, target_chembl_id):
    hits = activity.filter(
        molecule_chembl_id=molecule_chembl_id,
        target_chembl_id=target_chembl_id,
        standard_type="IC50",
        standard_units="nM",
        pchembl_value__isnull=False,
    ).only(["standard_value"])
    values = [float(h["standard_value"]) for h in hits if h["standard_value"] is not None]
    return min(values) if values else None


def get_target_activities(target_chembl_id, limit=2000):
    hits = activity.filter(
        target_chembl_id=target_chembl_id,
        standard_type="IC50",
        standard_units="nM",
        pchembl_value__isnull=False,
    ).only(["standard_value"])[:limit]
    return [float(h["standard_value"]) for h in hits if h["standard_value"] is not None]


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
    preds = model.predict(smiles=smiles)
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
    return list(hits)


def get_similar_smiles(similar_compounds):
    result = []
    for c in similar_compounds:
        raw = molecule.get(c["molecule_chembl_id"])["molecule_structures"]["canonical_smiles"]
        result.append((desalt(raw), c["molecule_chembl_id"]))
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
    similar = get_similar_compounds(smiles, threshold=70)
    pairs = get_similar_smiles(similar)
    deduped = dedup_against_query(pairs, smiles)
    all_compounds = [(smiles, chembl_id)] + deduped

    with open("compounds.smi", "w") as f:
        for smi, cid in all_compounds:
            f.write(f"{smi}\t{cid}\n")
    subprocess.run(["mmpdb", "fragment", "compounds.smi", "-o", "fragments.fragdb"], check=True)
    subprocess.run(["mmpdb", "index", "fragments.fragdb", "-o", "mmp.db"], check=True)

    batch_preds = model.predict(smiles=[s for s, _ in all_compounds])
    with open("props.txt", "w") as f:
        f.write(f"ID\t{property_name}\n")
        for (s, cid), val in zip(all_compounds, batch_preds[property_name]):
            f.write(f"{cid}\t{val}\n")
    subprocess.run(["mmpdb", "loadprops", "-p", "props.txt", "mmp.db"], check=True)


# ── Agent 3: 구조 변형 제안 ───────────────────────────────────
def run_agent3(smiles, property_name):
    subprocess.run(
        ["mmpdb", "transform", "--smiles", smiles, "mmp.db",
         "--property", property_name, "-o", "transform_results.csv"],
        check=True,
    )
    df = pd.read_csv("transform_results.csv", sep="\t")
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
    diag = run_agent0_1(drug_name)
    if diag is None:
        return {"status": "타겟_판단불가"}
    if diag["property_name"] is None:
        return {"status": "독성_문제_없음", "진단": diag}

    run_agent2(diag["smiles"], diag["chembl_id"], diag["property_name"])
    all_variants = run_agent3(diag["smiles"], diag["property_name"])

    checked = 0
    for n in top_n_steps:
        new_batch = all_variants.iloc[checked:n]
        checked = n
        if len(new_batch) == 0:
            continue
        ranked = run_agent4_5(diag["smiles"], new_batch, diag["property_name"])
        if len(ranked) > 0:
            return {
                "status": "완료",
                "타겟문제": diag["property_name"],
                "랭킹": ranked,
                "확인한_후보수": n,
            }

    return {"status": "통과_후보_없음", "타겟문제": diag["property_name"], "확인한_후보수": checked}


def run_batch(drug_names):
    # 여러 화합물 순차 실행, 하나 에러 나도 중단 없이 기록만 남기고 계속
    results = {}
    for name in drug_names:
        try:
            results[name] = run_full_pipeline(name)
        except Exception as e:
            results[name] = {"status": "에러", "message": str(e)}
    return results


if __name__ == "__main__":
    # 검증됐던 화합물 3개로 기본 동작 확인
    batch_results = run_batch(["Cisapride", "Troglitazone", "Bromfenac"])
    for name, r in batch_results.items():
        print(name, "->", r.get("status"), "| 문제축:", r.get("타겟문제"))
