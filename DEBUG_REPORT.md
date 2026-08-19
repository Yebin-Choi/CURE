# CURE 파이프라인 — 배치 테스트/디버깅 세션 리포트

날짜: 2026-08-19
브랜치: `claude/cure-pipeline-batch-debug-zhsi1h`
대상 파일: `cure_pipeline.py`

## 1. 요청 사항

1. 이 실행 환경에서 ChEMBL(`www.ebi.ac.uk`), Deep-PK(`biosig.lab.uq.edu.au`) 접속 확인
2. 검증된 3개(Cisapride, Troglitazone, Bromfenac) 외 7개 화합물로 `run_batch()` 확장 실행
   - hERG: Terfenadine, Astemizole
   - DILI: Trovafloxacin
   - 안전 대조군: Metformin, Amoxicillin, Loratadine, Ibuprofen
3. 에러 발생 시 화합물·함수·원인을 로그로 남기고 수정 후 재실행
4. 화합물별 결과 표 (상태 / 문제축 / 통과 후보 수 / 최상위 후보 SMILES·SAScore)
5. Agent 0~5 역할 구분 유지하며 안정성 개선

## 2. 네트워크 확인 결과 — 실험 실행 불가 확정

이 세션의 아웃바운드 프록시(`/root/.ccr/README.md` 기반 정책 프록시)에서 두 도메인 모두
**403(조직 정책 차단)**으로 막혀 있음을 확인했다. 세션 시작 시점과 세션 종료 시점 두 번
재확인했고 결과는 동일했다.

```
$ curl https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json
curl: (56) CONNECT tunnel failed, response 403

$ curl https://biosig.lab.uq.edu.au/deeppk/
curl: (56) CONNECT tunnel failed, response 403
```

프록시 상태 엔드포인트(`$HTTPS_PROXY/__agentproxy/status`)의 `recentRelayFailures`에도
두 호스트에 대한 `connect_rejected` (`gateway answered 403 to CONNECT (policy denial or
upstream failure)`) 기록이 남아 있다. 이는 일시적 장애가 아니라 **이 세션이 속한 워크스페이스의
아웃바운드 egress 정책 자체가 이 두 도메인을 차단**하는 것으로 판단된다 (PyPI 등 화이트리스트에
없는 임의 외부 호스트는 기본 차단인 것으로 보임).

**결론**: 요청 #2(배치 확장 실행), #3(실패 화합물 재실행), #4(결과 표)는 **이 실행 환경에서는
수행 불가**. 이 환경의 egress 정책에 두 도메인을 허용해주거나, ChEMBL/Deep-PK에 접근 가능한
다른 환경(Colab 등)에서 아래 "실행 방법" 절의 스크립트를 그대로 돌리면 된다. 우회 시도는 하지
않았다 (요청 원칙에 따름).

## 3. 코드에 적용한 변경 (`cure_pipeline.py`, 커밋 `02ef057`)

실제 배치를 못 돌리는 대신, 배경 설명에서 언급된 "실제로 겪었던" 실패 패턴들에 대한 방어 코드를
추가했다. Agent 0~5의 역할 구분과 함수 시그니처의 의미는 바꾸지 않았다 (Agent 3만 `chembl_id`
파라미터가 하나 늘었는데, 이는 화합물별 mmpdb 파일 경로를 찾기 위한 것이고 로직 자체는 동일).

| 변경 | 내용 | 이유 |
|---|---|---|
| `_with_retry()` | ChEMBL/ADMET-AI/Deep-PK 호출을 2/4/8/16초 지수 백오프로 최대 5회 재시도. 실패마다 어느 화합물·어느 호출인지 로그 출력, 최종 실패 시 그 정보를 담은 `RuntimeError` 발생 | 배경 설명: "500/404/타임아웃이 실제로 자주 발생" + "에러를 잡아서 어느 화합물의 어느 단계에서 실패했는지 로그로 남겨줄 것" |
| `_work_path()` | mmpdb 산출물(`compounds.smi`, `fragments.fragdb`, `mmp.db`, `props.txt`, `transform_results.csv`)을 `mmpdb_work/<chembl_id>/` 디렉터리로 분리 | 배경 설명: "화합물마다 덮어써짐... 배치로 여러 화합물을 돌릴 때... 파일명에 화합물 식별자를 넣는 방식으로 분리" |
| `_run_mmpdb()` | mmpdb subprocess 실패 시 stderr를 그대로 출력 후 raise | 기존엔 `CalledProcessError`만 던져서 mmpdb 내부 에러 메시지가 안 보였음 |
| `deeppk_predict()` | 응답에 예상 키가 없으면 **실제 응답 구조**를 담아 `KeyError` 발생 | 배경 설명: "원인 추정만 하지 말고 실제 응답/데이터 구조를 먼저 찍어서 확인" — Deep-PK 응답 스펙이 바뀌거나 항목명이 틀렸을 때 바로 알 수 있게 |
| `get_similar_smiles()` 버그 수정 | 재시도 람다가 루프 변수 `cid`를 늦은 바인딩(late binding)으로 캡처해 마지막 화합물 ID로 전부 덮어써지던 버그 수정 (`lambda cid=cid: ...`) | 재시도 로직 추가 과정에서 발견 |
| `run_full_pipeline()` | Agent 0-1 / 2 / 3 / 4-5 각 단계를 try/except로 감싸 `[화합물명] AgentX 단계 실패: ...` 형태로 재발생 | `run_batch()` 로그만 보고 어느 화합물이 어느 Agent에서 죽었는지 바로 알 수 있게 |
| `run_batch()` | 화합물 시작 시 `=== {name} 실행 시작 ===` 출력, 에러 시 즉시 `print` (기존엔 조용히 dict에만 저장) | 배치 실행 중 실시간 진행 상황 확인 |
| `__main__` 블록 | 배치 목록에 요청받은 7개 화합물 추가 | 네트워크 열린 환경에서 바로 요청 #2를 실행할 수 있도록 |

### 3.1 검증 방법 (오프라인)

`admet_ai`, `chembl_webresource_client`를 목(mock)으로 치환해 실제 네트워크 없이 순수 로직만
단위 테스트했다 (`rdkit`, `pandas`, `numpy`, `requests`는 실제 패키지 사용).

- `_with_retry`: N번 실패 후 성공 시 정상 반환 / 전량 실패 시 label이 포함된 `RuntimeError` 확인
- `_work_path`: 화합물별로 다른 경로 생성 + 디렉터리 실제 생성 확인
- `get_similar_smiles`: 클로저 버그 수정 후 cid가 순서대로 정확히 매칭되는지 확인
- `deeppk_predict`: 응답에 키가 없을 때 실제 응답 구조가 에러 메시지에 노출되는지 확인
- `poll_deeppk`: `running`(dict) → 완료(이중 인코딩 JSON 문자열) 흐름이 정상 파싱되는지 확인

6개 테스트 모두 통과. `python -m py_compile cure_pipeline.py`로 전체 문법도 확인했다.

**한계**: mmpdb CLI, ADMET-AI 실제 모델, 실제 ChEMBL/Deep-PK 응답 스펙은 이 환경에서 실행할 수
없어 오프라인 테스트로 검증하지 못했다. 특히 `DEEPPK_KEY_MAP`의 키 이름(`"[Toxicity/Liver
Injury I (DILI)] Predictions"` 등)이 실제 API 응답과 일치하는지는 네트워크가 열린 환경에서
Terfenadine/Astemizole(hERG)이나 Trovafloxacin(DILI)으로 첫 실행할 때 확인이 필요하다 — 여기서
막혔다면 이제는 `deeppk_predict`가 실제 응답 구조를 그대로 에러에 찍어주므로 원인 파악이 빠를
것이다.

## 4. 실행 방법 (네트워크 열린 환경)

```bash
pip install chembl_webresource_client admet-ai mmpdb rdkit requests pandas numpy
python cure_pipeline.py
```

`__main__` 블록이 아래 10개 화합물로 `run_batch()`를 실행하고 화합물별 상태/문제축을 출력한다.
전체 랭킹 데이터(`SMILES`, `SAScore`, 개선도 등)가 필요하면 `run_batch(...)`의 반환값을 그대로
받아서 각 화합물의 `["랭킹"]` DataFrame을 확인하면 된다 (요청하신 "상태 / 문제축 / 통과 후보 수 /
최상위 후보 SMILES·SAScore" 표는 이 DataFrame에서 바로 뽑을 수 있음: `df.iloc[0][["SMILES",
"SAScore"]]` + `len(df)`).

```python
from cure_pipeline import run_batch

names = [
    "Cisapride", "Troglitazone", "Bromfenac",
    "Terfenadine", "Astemizole", "Trovafloxacin",
    "Metformin", "Amoxicillin", "Loratadine", "Ibuprofen",
]
results = run_batch(names)

rows = []
for name, r in results.items():
    status = r.get("status")
    axis = r.get("타겟문제")
    if status == "완료":
        top = r["랭킹"].iloc[0]
        rows.append([name, status, axis, len(r["랭킹"]), top["SMILES"], top["SAScore"]])
    else:
        rows.append([name, status, axis, 0, None, None])

import pandas as pd
print(pd.DataFrame(rows, columns=["화합물", "상태", "문제축", "통과후보수", "최상위SMILES", "SAScore"]))
```

## 5. 이 세션에서 겪은 별도 이슈 (참고용)

작업과 별개로, 이 저장소에 대한 **Claude의 GitHub push 권한**이 세션 시작 시점에 없었다
(`403 Resource not accessible by integration`). claude.ai 커넥터의 GitHub 토글(OAuth 인증)만으로는
부족했고, 실제로는 `https://github.com/apps/claude`에서 **GitHub App을 이 저장소에 직접
설치**(Repository access에 `CURE` 추가, Contents: Read and write 권한 포함)해야 했다. 이후
정상적으로 브랜치 push가 됐다. 같은 문제를 겪는다면 이 경로부터 확인하는 게 빠르다.

## 6. 다음 단계 제안

1. 네트워크가 열린 환경(Colab, 로컬 등)에서 위 "실행 방법"으로 10개 화합물 배치 실행
2. 에러가 나면 `_with_retry`/`_run_mmpdb`/`deeppk_predict`가 남기는 로그(화합물명 + 실제 응답/에러
   구조 포함)를 보고 원인 파악 → 원인만 잡히면 그 화합물부터 재실행
3. 특히 hERG/DILI 케이스(Terfenadine, Astemizole, Trovafloxacin)는 `DEEPPK_KEY_MAP`을 처음
   타는 경로라 Deep-PK 응답 키 이름이 실제와 맞는지부터 확인 권장
4. 결과 표는 위 스니펫으로 바로 생성 가능
