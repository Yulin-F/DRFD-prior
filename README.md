# DRFD-prior

This repository contains the implementation of the paper: Distributionally Robust Fault Detection Trade-off Design with Prior Fault Information ([Arxiv](arxiv.org/abs/2412.20237)).

## Problem Formulation

We presented a new DRFD trade-off design scheme in this work by utilizing prior fault information. The key contribution includes a novel distributional robustness metric of detecting a known fault and a new relaxed distributionally robust chance constraint that ensures robust detectability. Then, a new DRFD design problem of fault detection under unknown probability distributions is proposed, and this offers a flexible balance between the robustness of detecting known critical faults and the overall detectability against all possible faults. 
$$
\max_{P,\eta>
0} &~ \rho\left(P\right) +  {\gamma}{\eta}  \\
{\rm s.t.} &~ \sup_{\mathbb{P}_{\xi} \in \mathcal{D}_{\rm W}(\theta, N)}
\mathbb{P}_\xi \left\{\lVert P \xi\rVert^2> 1 \right\} \leq \alpha\\
&~\sup_{\mathbb{P}_{\xi}  \in \mathcal{M}(\Xi)} \mathbb{P}_\xi \left\{\lVert  P \left( \xi + V \bar{f} \right) \rVert^2 \leq 1\right\}  -  {\frac{1}{\eta}} d_{\rm W}(\mathbb{P}_{\xi},\hat{\mathbb{P}}_{N})
  \leq 1 -\beta  ,\\
$$

## Get Start

1- Download [MOSEK](https://www.mosek.com) and [YALMIP]([YALMIP](https://yalmip.github.io/)).
2- Add MOSEK and YALMIP folder and subfolders to the Matlab path. 
3- Add the repo folder and subfolders to the Matlab path.  

## Citation

If you found this repository useful, please consider citing:

```tex
@article{feng2024distributionally,
  title={Distributionally Robust Fault Detection Trade-off Design with Prior Fault Information},
  author={Feng, Yulin and Jin, Hailang and Ding, Steven X and Ye, Hao and Shang, Chao},
  journal={arXiv preprint arXiv:2412.20237},
  year={2024}
}
```

## Contact us

To contact us about DRFD-prior, email either [Yulin Feng](mailto:fyl23@mails.tsinghua.edu.cn?Subject=DRFD-prior).





