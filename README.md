# S2S Dynamics Based Bipedal Walking

Slowly uploading the implementations on the S2S Dynamics based framework for controlling bipedal robotic walking. 

 ### cpp folder now contains basic examples of using the S2S based framework with 
 - a H-LIP class 
 - walking outputs using H-LIP approximation 
 - a QP based tracking controller
 
 TODO: A full repo of using it to control Cassie in Mujoco will be uploaded soon. 
 
 ### Matlab folder mainly contains a LIPapprox class (still un-organzied yet since it was research code) with controllers 
 - Deadbeat, LQR, CLF, etc
 - MPC, SLS
 - Adaptive Controller
 
 TODO: I will try to organize it in a better way. These linear controllers are relatively straight forward, only some of which have been appeared in previous publications.

---
Most walking videos can be seen from the [youtube channel](https://www.youtube.com/channel/UC__Fnw5l_TQCIBIIrfKXBtA).

Results of the S2S based framework: 
- 3D Periodic Walking: [Simulation](https://www.youtube.com/watch?v=-_QmNNBPfdg), [Experiment](https://www.youtube.com/watch?v=9DtRkHP_tQU)
- [**Unstable Terrain**](https://www.youtube.com/watch?v=DOS-xBs4Kdw)
- [Versatile Walking](https://www.youtube.com/watch?v=aFLNkHYTDaw)
- Rough Terrain: [aSLIP](https://www.youtube.com/watch?v=fUZu6y-Gu4g), [Cassie](https://www.youtube.com/watch?v=mHboC-vUZhM)   
- Push Rejection: [Sagittal](https://www.youtube.com/watch?v=MeaR__wgYyY), [Coronal](https://www.youtube.com/watch?v=_EqxuzywQWU)   
- [Global Path Tracking](https://www.youtube.com/watch?v=06efo-U1mrw) 

---
[Paper (Brief)](https://arxiv.org/pdf/2101.09588.pdf):
```
Xiaobin Xiong, and Aaron Ames. 
3D underactuated bipedal walking via h-lip based gait synthesis and stepping stabilization.
arXiv preprint arXiv:2101.09588 (2021).
```

[Thesis](https://thesis.library.caltech.edu/14230/1/Thesis__to_submit_0601.pdf)  
```
Xiong, Xiaobin
Reduced Order Model Inspired Robotic Bipedal Walking: A Step-to-step Dynamics Approximation based Approach
PHD Thesis, California Institute of Technology
```
