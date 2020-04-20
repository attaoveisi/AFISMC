# AFISMC
a new observer-based adaptive fuzzy integral sliding mode controller (AFISMC) is proposed based on the Lyapunov stability theorem. The plant under study is subjected to a square-integrable disturbance and is assumed to have mismatch uncertainties both in state- and input-matrices.  In addition, a norm-bounded time varying term is introduced to address the possible existence of un-modelled/nonlinear dynamics. Based on the classical sliding mode controller (SMC), the equivalent control effort is obtained to satisfy the sufficient requirement of SMC and then the control law is modified to guarantee the reachability of the system trajectory to the sliding manifold. The sliding surface is compensated based on the observed states in the form of linear matrix inequality (LMI). In order to relax the norm-bounded constrains on the control law and solve the chattering problem of SMC, a fuzzy logic (FL) inference mechanism is combined with the controller. An adaptive law is then introduced to tune the parameters of the fuzzy system on-line. Finally, by aiming at evaluating the validity of the controller and the robust performance of the closed-loop system, the proposed regulator is implemented on a real-time mechanical vibrating system.

Highlights:
A new observer-based adaptive fuzzy integral sliding mode controller (AFISMC) is proposed based on the Lyapunov stability theorem. Some important aspects of this paper are:

An appropriate switching function based on the observed states is designed to satisfy the classical sliding mode (SM) constraints.
An observer system is designed simultaneously based on the measured output and control signal by using Lyapunov stability theorem.
Matched modelling uncertainty representation for state and input matrices together with a norm-bounded un-modeled dynamics is considered in the design procedure to guarantee the stability of close-loop system.
By introduction of  disturbance rejection Lemma, the desired attenuation rate of the unknown norm-bounded time varying disturbance on measured output is achieved.
The conservatism of the control system due to the existence of various uncertain terms and disturbance is addressed by a fuzzy system. Introduction of the fuzzy system presented the following advantages:
The chattering problem of the SMC is addressee properly.
The reaching phase of the trajectory of the system on the SM manifold is guaranteed.
By addition of an adaptive intelligent algorithm base on online tuning of fuzzy parameters, the conservatism problem of the ISMC is relieved. 
Despite most of the cited researches, a real time application of the proposed control system is presented to investigate the different aspects of the performance of the closed-loop system such as:

Disturbance rejection quality in the nominal frequency range of the modelling is presented.
Robustness of the performance is investigated with respect to higher order un-modeled dynamics.
Sensitivity of the applied control effort is studied with respect to spillover effect.
