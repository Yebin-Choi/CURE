"""
오프라인 모의 실험 하네스.

ChEMBL/ADMET-AI/Deep-PK 호출부만 실제와 비슷한 합성 데이터로 대체하고,
RDKit(구조 처리)·mmpdb(실제 CLI 서브프로세스)·pandas 로직은 cure_pipeline.py의
진짜 코드를 그대로 실행한다. 목적은 네트워크 없이도 여러 화합물을 넣어봤을 때
파이프라인 로직 자체(파일 처리, mmpdb 연동, 게이트/랭킹 로직, 예외 처리)에
버그가 있는지 찾아내는 것.
"""
import hashlib
import os
import random
import sys
import types

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ── admet_ai / chembl_webresource_client를 실제 모듈처럼 흉내내는 모의 계층 ──
admet_ai_mod = types.ModuleType("admet_ai")


TOX_ENDPOINTS = [
    "AMES", "DILI", "hERG", "ClinTox", "Carcinogens_Lagunin", "Skin_Reaction",
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
    "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
]


def _seeded_score(smiles, endpoint, base=None):
    h = hashlib.md5(f"{smiles}|{endpoint}".encode()).hexdigest()
    v = int(h[:8], 16) / 0xFFFFFFFF
    if base is not None:
        v = min(1.0, max(0.0, base + (v - 0.5) * 0.15))
    return round(v, 3)


# 시나리오별로 특정 화합물엔 특정 축이 "고위험"이 되도록 base를 강제한다
FORCED_PRIMARY = {
    "Terfenadine": ("hERG", 0.9),
    "Astemizole": ("hERG", 0.88),
    "Trovafloxacin": ("DILI", 0.85),
    "Cisapride": ("hERG", 0.82),
    "Troglitazone": ("DILI", 0.81),
    "Bromfenac": ("AMES", 0.78),  # DEEPPK_KEY_MAP에 없는 축 경로 테스트
    "HyphenTestDrug": ("NR-AR-LBD", 0.75),  # 하이픈 포함 property명이 mmpdb를 통과하는지 테스트
}
SAFE_CONTROLS = {"Metformin", "Amoxicillin", "Loratadine", "Ibuprofen"}


class FakeADMETModel:
    def __init__(self, *a, **k):
        pass

    def predict(self, smiles):
        if isinstance(smiles, str):
            return self._predict_one(smiles, drug_hint=getattr(self, "_current_drug", None))
        # 배치: {endpoint: [값들]} 형태, 입력 순서와 정렬 일치해야 함
        rows = [self._predict_one(s) for s in smiles]
        return {ep: [r[ep] for r in rows] for ep in TOX_ENDPOINTS}

    def _predict_one(self, smiles, drug_hint=None):
        forced_ep, forced_base = (None, None)
        if drug_hint and drug_hint in FORCED_PRIMARY:
            forced_ep, forced_base = FORCED_PRIMARY[drug_hint]
        row = {}
        for ep in TOX_ENDPOINTS:
            if drug_hint in SAFE_CONTROLS:
                row[ep] = _seeded_score(smiles, ep, base=0.15)
            elif ep == forced_ep:
                row[ep] = forced_base
            else:
                row[ep] = _seeded_score(smiles, ep, base=0.3)
        return row


admet_ai_mod.ADMETModel = FakeADMETModel
sys.modules["admet_ai"] = admet_ai_mod

chembl_pkg = types.ModuleType("chembl_webresource_client")
chembl_new_client_mod = types.ModuleType("chembl_webresource_client.new_client")


class FakeQS(list):
    def only(self, fields):
        return self

    def filter(self, **kw):
        return self


class FakeEndpoint:
    def __init__(self, kind, registry):
        self.kind = kind
        self.registry = registry

    def search(self, name):
        entry = self.registry.get(name)
        if entry is None or entry.get("no_chembl_id"):
            return FakeQS([])
        return FakeQS([{"molecule_chembl_id": entry["chembl_id"]}])

    def filter(self, **kw):
        if self.kind == "mechanism":
            cid = kw.get("molecule_chembl_id")
            entry = self._by_cid(cid)
            if entry is None or entry.get("no_target"):
                return FakeQS([])
            return FakeQS([{"target_chembl_id": entry["target_id"]}])
        if self.kind == "activity":
            cid = kw.get("molecule_chembl_id")
            tid = kw.get("target_chembl_id")
            if cid is not None:
                entry = self._by_cid(cid)
                if entry is None or entry.get("insufficient_activity"):
                    return FakeQS([])
                return FakeQS([{"standard_value": str(entry.get("own_ic50", 50.0))}])
            # target 전체 population
            entry = self._by_target(tid)
            if entry is None:
                return FakeQS([])
            pop = entry.get("population_ic50", [])
            return FakeQS([{"standard_value": str(v)} for v in pop])
        if self.kind == "similarity":
            smi = kw.get("smiles")
            entry = self._by_smiles(smi)
            if entry is None:
                return FakeQS([])
            return FakeQS(
                [{"molecule_chembl_id": cid, "similarity": 80} for cid, _ in entry.get("analogs", [])]
            )
        return FakeQS([])

    def get(self, chembl_id):
        own_entry = self._by_cid(chembl_id)
        if own_entry is not None:
            return {"molecule_structures": {"canonical_smiles": own_entry["smiles"]}}
        return {"molecule_structures": {"canonical_smiles": self._analog_smiles(chembl_id)}}

    def _by_cid(self, cid):
        for e in self.registry.values():
            if isinstance(e, dict) and e.get("chembl_id") == cid:
                return e
        return None

    def _by_target(self, tid):
        for e in self.registry.values():
            if isinstance(e, dict) and e.get("target_id") == tid:
                return e
        return None

    def _by_smiles(self, smi):
        for e in self.registry.values():
            if isinstance(e, dict) and e.get("smiles") == smi:
                return e
        return None

    def _by_analog_cid(self, cid):
        for e in self.registry.values():
            if not isinstance(e, dict):
                continue
            for acid, _ in e.get("analogs", []):
                if acid == cid:
                    return e
        return None

    def _analog_smiles(self, cid):
        for e in self.registry.values():
            if not isinstance(e, dict):
                continue
            for acid, asmi in e.get("analogs", []):
                if acid == cid:
                    return asmi
        raise KeyError(cid)


class FakeNewClient:
    def __init__(self, registry):
        self.molecule = FakeEndpoint("molecule", registry)
        self.mechanism = FakeEndpoint("mechanism", registry)
        self.target = FakeEndpoint("target", registry)
        self.activity = FakeEndpoint("activity", registry)
        self.similarity = FakeEndpoint("similarity", registry)


REGISTRY = {}
chembl_new_client_mod.new_client = FakeNewClient(REGISTRY)
chembl_pkg.new_client = chembl_new_client_mod
sys.modules["chembl_webresource_client"] = chembl_pkg
sys.modules["chembl_webresource_client.new_client"] = chembl_new_client_mod

import cure_pipeline as cp  # noqa: E402

cp.time.sleep = lambda s: None


# ── analog(유사 화합물) 생성: 단순 치환으로 mmpdb가 실제 매칭쌍을 찾을 수 있게 함 ──
def make_analogs(smiles, n=12):
    from rdkit import Chem

    edits = [
        ("Cl", "F"), ("F", "Cl"), ("Br", "Cl"), ("Br", "F"),
        ("OC", "OCC"), ("OCC", "OCCC"), ("CC(C)C", "CC(C)CC"),
        ("c1ccccc1", "c1ccc(F)cc1"), ("c1ccccc1", "c1ccc(Cl)cc1"),
        ("C(=O)O", "C(=O)OC"), ("N", "NC"), ("N2CC", "N2CCC"),
        ("CC1", "CCC1"), ("c1ccc(", "c1cc(F)c("), ("cc1", "c(C)c1"),
    ]
    canon0 = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    seen = {canon0}
    out = []
    random.seed(hash(smiles) & 0xFFFF)
    for old, new in edits:
        if len(out) >= n:
            break
        idx = 0
        while old in smiles[idx:] and len(out) < n:
            pos = smiles.index(old, idx)
            cand = smiles[:pos] + new + smiles[pos + len(old):]
            idx = pos + len(old)
            mol = Chem.MolFromSmiles(cand)
            if mol is None:
                continue
            canon = Chem.MolToSmiles(mol)
            if canon not in seen:
                seen.add(canon)
                out.append(canon)
    return out


# ── 테스트 화합물 정의 (실제로 알려진 SMILES, 정밀도보다 코드 검증 목적) ──
RAW_COMPOUNDS = {
    "Cisapride": "COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC",
    "Troglitazone": "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
    "Bromfenac": "Nc1ccc(C(=O)c2ccc(Br)cc2)c(CC(=O)O)c1",
    "Terfenadine": "CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",
    "Astemizole": "COc1ccc(CCN2CCC(Nc3nc4ccc(F)cc4n3Cc3ccc(F)cc3)CC2)cc1",
    "Trovafloxacin": "OC(=O)c1cn(C2CC2)c2cc(N3CC4CC3CN4)c(F)cc2c1=O",
    "Metformin": "CN(C)C(=N)NC(=N)N",
    "Amoxicillin": "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
    "Loratadine": "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "HyphenTestDrug": "COc1ccc(C(=O)Nc2ccc(Cl)cc2)cc1N",
    # 엣지 케이스
    "존재하지않는약": None,           # ChEMBL ID 조회 실패 -> None
    "TargetlessDrug": "CCO",           # 타겟 미확정 케이스용 (에탄올로 단순화)
}

EDGE_NO_CHEMBL = {"존재하지않는약"}
EDGE_NO_TARGET = {"TargetlessDrug"}

def _canon(smiles):
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


for i, (name, smi) in enumerate(RAW_COMPOUNDS.items()):
    if name in EDGE_NO_CHEMBL:
        REGISTRY[name] = {"no_chembl_id": True}
        continue
    smi = _canon(smi)  # get_smiles()가 desalt()로 재정규화하므로 등록 시점에 맞춰둔다
    cid = f"CHEMBLFAKE{i:03d}"
    tid = f"TGTFAKE{i:03d}"
    entry = {
        "chembl_id": cid,
        "target_id": tid,
        "smiles": smi,
        "own_ic50": 40.0 + i * 5,
        "population_ic50": [10, 20, 30, 40, 60, 80, 120, 200, 300][: 5 + (i % 4)],
    }
    if name in EDGE_NO_TARGET:
        entry["no_target"] = True
    analogs = make_analogs(smi, n=12)
    entry["analogs"] = [(f"{cid}_A{j}", a) for j, a in enumerate(analogs)]
    REGISTRY[name] = entry

# FakeADMETModel이 "지금 어느 화합물을 보고 있는지" 알 수 있도록 run_agent0_1 호출 전에 힌트 주입
_orig_run_agent0_1 = cp.run_agent0_1


def patched_run_agent0_1(drug_name):
    cp.model._current_drug = drug_name
    return _orig_run_agent0_1(drug_name)


cp.run_agent0_1 = patched_run_agent0_1

# Deep-PK는 HTTP 형태까지 흉내내지 않고 결과 함수 자체를 모의 처리 (HTTP 파싱은 이미 별도 단위테스트로 검증됨)
_deeppk_call_count = {"n": 0}


def fake_deeppk_predict(smiles, key):
    _deeppk_call_count["n"] += 1
    # smiles 해시 기반으로 절반 정도는 Safe, 절반은 Toxic이 나오게
    h = int(hashlib.md5(smiles.encode()).hexdigest()[:4], 16)
    return "Safe" if h % 2 == 0 else "Toxic"


cp.deeppk_predict = fake_deeppk_predict

if __name__ == "__main__":
    names = list(RAW_COMPOUNDS.keys())
    print(f"=== 모의 배치 실행: {len(names)}개 화합물 ===\n")
    results = cp.run_batch(names)

    print("\n=== 결과 표 ===")
    rows = []
    for name, r in results.items():
        status = r.get("status")
        axis = r.get("타겟문제")
        checked = r.get("확인한_후보수", 0)
        if status == "완료":
            top = r["랭킹"].iloc[0]
            rows.append((name, status, axis, len(r["랭킹"]), checked, top["SMILES"], round(top["SAScore"], 2)))
        else:
            rows.append((name, status, axis, 0, checked, "-", "-"))
    width = max(len(r[0]) for r in rows)
    for row in rows:
        print(f"{row[0]:<{width}} | {row[1]:<12} | {str(row[2]):<10} | 통과 {row[3]:>3} | 확인 {row[4]:>3} | {row[5][:40]:<40} | SA {row[6]}")

    print(f"\nDeep-PK 모의 호출 횟수: {_deeppk_call_count['n']}")
