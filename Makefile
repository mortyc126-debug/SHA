# Makefile for SHA-256 GPU differential cryptanalysis experiments
#
# Targets:
#   all             - build all experiments
#   birthday        - birthday_search_17 (main A1 attack)
#   tower           - p_adic_tower (H1 hypothesis)
#   clean           - remove binaries
#
# GPU architecture flags:
#   sm_75 = Turing  (RTX 2080, T4)
#   sm_80 = Ampere  (A100, RTX 3090)
#   sm_86 = Ampere  (RTX 3080, A10)
#   sm_89 = Ada     (RTX 4090, L40)
#   sm_90 = Hopper  (H100)
#
# Auto-detect with: nvcc --run-dir . -arch=native
# Or specify: make ARCH=sm_89

NVCC     = nvcc
ARCH    ?= sm_80
NVCCFLAGS = -O3 -arch=$(ARCH) --use_fast_math -lineinfo \
            -Xcompiler "-O3 -march=native"

.PHONY: all birthday tower clean

all: birthday tower

birthday: birthday_search_17
p_adic: p_adic_tower

birthday_search_17: birthday_search_17.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<
	@echo "Built: $@"

p_adic_tower: p_adic_tower.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<
	@echo "Built: $@"

clean:
	rm -f birthday_search_17 p_adic_tower *.o

# Quick test: verify known pair П-15
test: birthday_search_17
	./birthday_search_17 e82222c7 1 1 test_pair.txt
	@echo "If pair W0=e82222c7,w1=516cfb41 appears -> PASS"

# Detect GPU architecture automatically
detect_arch:
	@python3 -c "import subprocess,re; \
	    r=subprocess.run(['nvcc','--version'],capture_output=True,text=True); \
	    print('nvcc:', r.stdout.split('release')[-1].split(',')[0].strip())"
	@nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
	    | head -1 | awk '{print "GPU arch: sm_" int($$1*10)}'
