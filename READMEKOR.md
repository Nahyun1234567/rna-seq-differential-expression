# RNA-seq 차등 발현 분석 (Differential Expression Analysis)

이 프로젝트는 Python과 통계적 방법을 활용하여 질병과 연관된 유전자를 식별하는  
엔드 투 엔드(end-to-end) RNA-seq 차등 발현 분석 파이프라인을 구현하는 것을 목표로 합니다.

---

## 프로젝트 개요 (Overview)

이 프로젝트는 장수와 관련된 특정 유전자 변이가 간암 세포의 유전자 발현 패턴에 어떤 영향을 미치는지를  
RNA-seq 데이터 분석을 통해 조사합니다.

최근 연구에 따르면, 100세 이상 장수한 사람들에게서 발견된 SIRT6 유전자 변이(N308K/A313S)는  
DNA 복구 능력을 향상시키고 암을 억제하는 특성과 연관되어 있습니다.  
하지만 이러한 변이가 실제 암 세포에서 어떤 분자적 메커니즘을 통해 작용하는지는 아직 명확히 밝혀지지 않았습니다.

본 프로젝트에서는  
야생형(wild-type) SIRT6 또는 장수인에서 발견된 SIRT6 변이를 발현하는 간암 세포의 RNA-seq 데이터를 분석하고,  
이를 대조군(control)과 비교합니다.

---

## 데이터셋 (Dataset)

- 출처: GEO (Gene Expression Omnibus)
- 데이터 유형: RNA-seq count matrix
- 샘플 구성: 질병군 vs 대조군
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314110

---

### 데이터셋 상세 설명 (Dataset Description)

- 데이터 출처: NCBI Gene Expression Omnibus (GEO), GSE314110  
- 생물 종: 인간 (Homo sapiens)  
- 실험 유형: RNA-seq (차세대 염기서열 분석 기반 유전자 발현 분석)

#### 실험 설계 (Experimental Design)

본 데이터셋은 두 가지 인간 간암 세포주를 포함합니다:
- HepG2
- Huh-7

각 세포주는 다음 세 가지 조건으로 처리되었습니다:
1. 대조군 (Control): Luciferase만 발현
2. SIRT6 WT: 일반적인(야생형) 인간 SIRT6 유전자 발현
3. SIRT6 변이형 (Variant): 장수인에서 발견된 SIRT6 변이(N308K/A313S) 발현

각 조건마다 3개의 생물학적 반복 샘플이 포함되어 있으며,  
총 18개의 RNA-seq 샘플로 구성되어 있습니다.

---

## 분석 방법 (Methods)

- 원시 count 데이터의 로그 정규화(Log normalization)
- t-test를 이용한 차등 발현 유전자 분석
- Benjamini–Hochberg 방법을 이용한 다중 가설 검정 보정(FDR)
- Volcano plot 및 Heatmap을 활용한 시각화

---

## 프로젝트 목표 (Goal)

- 대조군과 SIRT6 처리 간암 세포 간의 차등 발현 유전자(Differentially Expressed Genes, DEGs)를 식별
- 장수인에서 발견된 SIRT6 변이가  
  암의 공격성 및 침습성과 관련된 유전자 발현을 억제하는지를 분석

---

## 저장소 구조 (Repository Structure)

- `data/` : 원본 및 전처리된 데이터
- `notebooks/` : 단계별 분석 노트북
- `src/` : 재사용 가능한 Python 모듈
- `results/` : 시각화 결과 및 출력 테이블
- `report/` : 분석 결과 요약

---

## 활용 기술 (Skills Demonstrated)

- 생물정보학 데이터 분석
- 통계적 가설 검정
- 데이터 시각화
- 재현 가능한 연구 워크플로우 구축

---

### 한 줄 요약

> 장수 유전자 변이가 간암 세포의 유전자 발현을 어떻게 변화시키는지를 분석하는  
> RNA-seq 기반 Bioinformatics 프로젝트
