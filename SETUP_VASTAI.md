# Инструкция: запуск GPU-экспериментов на vast.ai

## Шаг 1 — Генерация SSH-ключа (PowerShell, Windows)

Открой PowerShell (Win+X → Windows PowerShell):

```powershell
# Генерация Ed25519 ключа (быстрый, рекомендуемый)
ssh-keygen -t ed25519 -C "vastai-sha256" -f "$env:USERPROFILE\.ssh\vastai_key"

# Когда спросит passphrase — можно оставить пустым (Enter дважды)

# Показать публичный ключ (он нужен для vast.ai)
Get-Content "$env:USERPROFILE\.ssh\vastai_key.pub"
```

Скопируй вывод (одна строка вида `ssh-ed25519 AAAA... vastai-sha256`).

---

## Шаг 2 — Добавление ключа на vast.ai

1. Войди на [vast.ai](https://vast.ai) → Account → SSH Keys
2. Нажми **"Add SSH Key"**
3. Вставь скопированный публичный ключ → Save

---

## Шаг 3 — Аренда GPU-инстанса

На странице **"Create"** на vast.ai выбери:

**Рекомендуемые параметры поиска:**
```
compute_cap >= 750          (GTX 1080 Ti, RTX 2080+, A100, ...)
cpu_arch in ['amd64']       (x86_64 — стандарт)
cuda_max_good >= 13         (CUDA 12+)
RAM >= 8 GB
```

**Рекомендуемые GPU (по соотношению цена/скорость):**
| GPU | SHA/s ожид. | Время 2^32 | Цена ~$/ч |
|-----|-------------|------------|-----------|
| A100 80GB | 110G | 0.04 сек | ~2.5 |
| RTX 4090  | 80G  | 0.05 сек | ~0.8 |
| RTX 3090  | 35G  | 0.12 сек | ~0.4 |
| A10       | 25G  | 0.17 сек | ~0.3 |
| T4        | 8G   | 0.5 сек  | ~0.1 |

**Шаблон запуска (On-Start Script):**

В поле "On-Start Script" (или `env`) вставь:
```bash
env >> /etc/environment
mkdir -p ${DATA_DIRECTORY:-/workspace}
cd ${DATA_DIRECTORY:-/workspace}
apt-get update -qq && apt-get install -y -qq git make 2>/dev/null
```

**Выбери Docker-образ:**
```
nvidia/cuda:12.3.2-devel-ubuntu22.04
```
(любой cuda devel образ с CUDA 12+)

---

## Шаг 4 — Подключение к инстансу

После запуска инстанса в vast.ai появится команда подключения вида:
```
ssh -p 12345 root@xxx.xxx.xxx.xxx
```

**В PowerShell:**
```powershell
# Подключиться (замени порт и адрес на реальные из vast.ai)
ssh -p PORT -i "$env:USERPROFILE\.ssh\vastai_key" root@IP_ADDRESS

# Если ключ не подхватывается автоматически:
ssh -p PORT -i "$env:USERPROFILE\.ssh\vastai_key" -o StrictHostKeyChecking=no root@IP_ADDRESS
```

---

## Шаг 5 — Загрузка кода на GPU-машину

**Вариант A: через git (рекомендуется)**
```bash
# На GPU-машине:
cd /workspace
git clone https://github.com/mortyc126-debug/SHA.git
cd SHA
git checkout claude/cuda-birthday-search-HTbNO
ls *.cu *.md Makefile
```

**Вариант B: через SCP (если git не работает)**
```powershell
# В PowerShell на Windows:
scp -P PORT -i "$env:USERPROFILE\.ssh\vastai_key" `
    birthday_search_17.cu p_adic_tower.cu Makefile run_experiments.sh `
    root@IP_ADDRESS:/workspace/
```

---

## Шаг 6 — Сборка

```bash
# На GPU-машине:
cd /workspace/SHA  # или /workspace если через SCP

# Определить архитектуру GPU автоматически
SM=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
echo "GPU arch: sm_$SM"

# Сборка
make all ARCH=sm_$SM

# Верификация (должна найти пару П-15: W0=e82222c7, w1=516cfb41)
make test
```

Ожидаемый вывод `make test`:
```
[Self-test] Verifying known pair П-15: ... OK
FOUND pair #1: w1=0x516cfb41  Da13=0x7711498a  DW16=0x84752d8e
```

---

## Шаг 7 — Запуск экспериментов

```bash
# Все эксперименты сразу (рекомендуется для первого запуска):
chmod +x run_experiments.sh
./run_experiments.sh 2>&1 | tee experiment_log.txt

# Или отдельно:

# --- A1: Birthday search (1 pair, быстро) ---
./birthday_search_17 e82222c7 1 1 test.csv
# Ожидается: найдёт пару П-15 за ~0.04 сек на A100

# --- A1: Collect 100 pairs for statistics (H2) ---
./birthday_search_17 0 1 100 pairs_100.csv
# 100 pairs: ~4 сек на A100, ~30 сек на RTX 3090

# --- H1: p-adic tower height test ---
./p_adic_tower 5000000 25 tower.txt
# 5M seeds, height up to 25: ~10 сек на A100

# --- Full 2^32 scan with known W0 (reproduces П-15 exactly) ---
./birthday_search_17 e82222c7 1 1 p15.csv
```

---

## Шаг 8 — Скачать результаты

```powershell
# В PowerShell на Windows (замените PORT и IP):
scp -P PORT -i "$env:USERPROFILE\.ssh\vastai_key" `
    "root@IP_ADDRESS:/workspace/SHA/results_*/*" `
    C:\Users\ИМЯ\Downloads\sha_results\
```

---

## Решение типичных проблем

**nvcc не найден:**
```bash
export PATH=/usr/local/cuda/bin:$PATH
which nvcc && nvcc --version
```

**Ошибка `sm_XY not supported`:**
```bash
# Проверить версию CUDA и поддерживаемые архитектуры:
nvcc --list-gpu-arch
# Собрать с явным указанием:
make ARCH=sm_86
```

**CUDA out of memory:**
```bash
# Уменьшить MAX_RESULTS в birthday_search_17.cu до 100000
# Или использовать меньший CHUNK (1ULL << 24 вместо 1ULL << 26)
```

**Медленная компиляция:**
```bash
# Добавить флаг для быстрой компиляции (меньше оптимизаций):
nvcc -O2 -arch=sm_$SM -o birthday_search_17 birthday_search_17.cu
```

---

## Ожидаемая производительность

| Эксперимент | A100 | RTX 4090 | RTX 3090 |
|-------------|------|----------|---------|
| 1 pair (2^32) | 0.04 сек | 0.05 сек | 0.12 сек |
| 100 pairs | ~4 сек | ~5 сек | ~12 сек |
| p_adic tower 1M seeds | ~2 сек | ~3 сек | ~7 сек |
| DW0 j-scan (32 values) | ~30 сек | ~40 сек | ~90 сек |

---

## Минимальная инструкция (быстрый старт)

```bash
# 1. Подключиться к GPU-машине
ssh -p PORT -i ~/.ssh/vastai_key root@IP

# 2. Клонировать и собрать
cd /workspace && git clone REPO_URL SHA && cd SHA
git checkout claude/cuda-birthday-search-HTbNO
SM=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
make ARCH=sm_$SM

# 3. Запустить
./birthday_search_17 e82222c7 1 1 test.csv   # должна найти П-15
./run_experiments.sh                           # все эксперименты
```
