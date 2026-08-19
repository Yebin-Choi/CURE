"""
mock_harness.py의 ChEMBL mock은 그대로 쓰고, ADMET-AI만 진짜 패키지로 교체해서
실제 화합물 구조에 대한 진짜 독성 예측값으로 전체 파이프라인을 검증한다.
(ChEMBL/Deep-PK는 여전히 네트워크 차단이라 mock 유지)

실행 전 `pip install admet-ai mmpdb rdkit pandas requests` 필요 (admet-ai는 완전
오프라인으로 동작 — 네트워크 불필요, 모델 가중치가 패키지에 내장돼 있음).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 1) 진짜 admet_ai를 그대로 두고, chembl_webresource_client만 목으로 교체
import types

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
        if entry is None:
            return FakeQS([])
        return FakeQS([{"molecule_chembl_id": entry["chembl_id"]}])

    def filter(self, **kw):
        if self.kind == "mechanism":
            cid = kw.get("molecule_chembl_id")
            entry = self._by_cid(cid)
            if entry is None:
                return FakeQS([])
            return FakeQS([{"target_chembl_id": entry["target_id"]}])
        if self.kind == "activity":
            cid = kw.get("molecule_chembl_id")
            tid = kw.get("target_chembl_id")
            if cid is not None:
                entry = self._by_cid(cid)
                if entry is None:
                    return FakeQS([])
                return FakeQS([{"standard_value": str(entry.get("own_ic50", 50.0))}])
            entry = self._by_target(tid)
            if entry is None:
                return FakeQS([])
            return FakeQS([{"standard_value": str(v)} for v in entry.get("population_ic50", [])])
        if self.kind == "similarity":
            entry = self._by_smiles(kw.get("smiles"))
            if entry is None:
                return FakeQS([])
            return FakeQS(
                [{"molecule_chembl_id": cid, "similarity": 80} for cid, _ in entry.get("analogs", [])]
            )
        return FakeQS([])

    def get(self, chembl_id):
        own = self._by_cid(chembl_id)
        if own is not None:
            return {"molecule_structures": {"canonical_smiles": own["smiles"]}}
        return {"molecule_structures": {"canonical_smiles": self._analog_smiles(chembl_id)}}

    def _by_cid(self, cid):
        for e in self.registry.values():
            if e.get("chembl_id") == cid:
                return e
        return None

    def _by_target(self, tid):
        for e in self.registry.values():
            if e.get("target_id") == tid:
                return e
        return None

    def _by_smiles(self, smi):
        for e in self.registry.values():
            if e.get("smiles") == smi:
                return e
        return None

    def _analog_smiles(self, cid):
        for e in self.registry.values():
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

import cure_pipeline as cp  # noqa: E402  (진짜 admet_ai를 그대로 import함)

cp.time.sleep = lambda s: None


def make_analogs(smiles, n=12):
    from rdkit import Chem

    edits = [
        ("Cl", "F"), ("F", "Cl"), ("Br", "Cl"), ("Br", "F"),
        ("OC", "OCC"), ("OCC", "OCCC"), ("CC(C)C", "CC(C)CC"),
        ("c1ccccc1", "c1ccc(F)cc1"), ("c1ccccc1", "c1ccc(Cl)cc1"),
        ("C(=O)O", "C(=O)OC"), ("N", "NC"), ("N2CC", "N2CCC"),
        ("CC1", "CCC1"), ("c1ccc(", "c1cc(F)c("), ("cc1", "c(C)c1"),
    ]
    from rdkit import Chem as _C
    canon0 = _C.MolToSmiles(_C.MolFromSmiles(smiles))
    seen = {canon0}
    out = []
    for old, new in edits:
        if len(out) >= n:
            break
        idx = 0
        while old in smiles[idx:] and len(out) < n:
            pos = smiles.index(old, idx)
            cand = smiles[:pos] + new + smiles[pos + len(old):]
            idx = pos + len(old)
            mol = _C.MolFromSmiles(cand)
            if mol is None:
                continue
            canon = _C.MolToSmiles(mol)
            if canon not in seen:
                seen.add(canon)
                out.append(canon)
    return out


def _canon(smiles):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


RAW = {
    "Cisapride": "COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC",
    "Terfenadine": "CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",
    "Trovafloxacin": "OC(=O)c1cn(C2CC2)c2cc(N3CC4CC3CN4)c(F)cc2c1=O",
    "Metformin": "CN(C)C(=N)NC(=N)N",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
}

for i, (name, smi) in enumerate(RAW.items()):
    smi = _canon(smi)
    cid = f"REALFAKE{i:03d}"
    tid = f"RTGT{i:03d}"
    analogs = make_analogs(smi, n=10)
    REGISTRY[name] = {
        "chembl_id": cid,
        "target_id": tid,
        "smiles": smi,
        "own_ic50": 40.0,
        "population_ic50": [10, 20, 30, 40, 60, 80, 120],
        "analogs": [(f"{cid}_A{j}", a) for j, a in enumerate(analogs)],
    }

import hashlib  # noqa: E402


def fake_deeppk_predict(smiles, key):
    h = int(hashlib.md5(smiles.encode()).hexdigest()[:4], 16)
    return "Safe" if h % 2 == 0 else "Toxic"


cp.deeppk_predict = fake_deeppk_predict  # Deep-PK도 여전히 네트워크 차단이라 mock 유지

if __name__ == "__main__":
    print("=== 진짜 ADMET-AI + mock ChEMBL/Deep-PK로 5개 화합물 실행 ===\n")
    for name in RAW:
        diag = cp.run_agent0_1(name)
        print(f"{name}: property_name(1차 문제축)={diag['property_name']}  "
              f"물성판정={diag['물성']['물성_판정']}  효능판정={diag['효능']['효능_판정']}")
    print()

    results = cp.run_batch(list(RAW.keys()))
    print("\n=== 최종 결과 ===")
    for name, r in results.items():
        status = r.get("status")
        axis = r.get("타겟문제")
        checked = r.get("확인한_후보수", 0)
        n_pass = len(r["랭킹"]) if status == "완료" else 0
        print(f"{name:<15} | {status:<12} | 축={str(axis):<10} | 확인={checked:>3} | 통과={n_pass}")
