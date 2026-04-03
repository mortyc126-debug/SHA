# Wang+Q∩T Hybrid: α-kernel Analysis

## Formula
α-kernel = max(0, 32*(R-8)) for R ≥ 9

## Data Table
| R | Ch linearized | Maj remaining | Filtered kernel | α-rank | α-kernel |
|:-:|:---:|:---:|:---:|:---:|:---:|
| 6 | 192 | 128 | 64 | 64 | 0 |
| 7 | 224 | 160 | 128 | 128 | 0 |
| 8 | 256 | 192 | 192 | 192 | 0 |
| **9** | 288 | 224 | 256 | 224 | **32** |
| **10** | 320 | 256 | 320 | 256 | **64** |
| **12** | 384 | 320 | 448 | 320 | **128** |
| **16** | 512 | 448 | 704 | 448 | **256** |

## Method
1. Wang chain: de[2..R]=0, gives e[r] values as constants
2. f[r]=e[r-1], g[r]=e[r-2] also known
3. Ch(e,f,g) products → constants (linearized)
4. Only Maj(a,b,c) products remain quadratic
5. α-system (quad eqs in kernel coordinates) has rank < dim at R≥9
