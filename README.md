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

Refs:

[1] J.Y. Hung, W.B. Gao, J.C. Hung, Variable structure control: A survey, IEEE Trans. Ind. Electron. 40 (1993) 2-22.

[2] Q.P. Ha, D.C. Rye, H.F. Durrant-Whyte, Robust sliding mode control with application, Int. J. Control. 72 (1999)
1078-1096.

[3] J. Ackermann, V.I. Utkin, Sliding mode control design based on Ackermann’s formula, IEEE T. Automat. Contr. 43(1998) 234-237.

[4] V.I. Utkin, J. Shi, Integral sliding mode in systems operating under uncertainty conditions, Proceedings of the 35th Conference on Decision and Control, Kobe, Japan, (1996) 4591-4596.

[5] C. Mnasri, M. Gasmi, LMI–based adaptive fuzzy integral sliding mode control of mismatched uncertain systems, Int. J. Appl. Math. Comput. Sci. 21 (2011) 605-615.

[6] W.J. Cao, J.X. Xu, Nonlinear integral-type sliding surface for both matched and unmatched uncertain systems, IEEE T. Automat. Contr. 49 (2004) 1355-1360.

[7] Q. Shaocheng, W. Yongji, Robust control of uncertain time delay system: A novel sliding mode control design via LMI, J. Syst. Eng. Electron. 17 (2006) 624-628.

[8] H. H. Choi, LMI-based sliding surface design for integral sliding mode control of mismatched uncertain systems, IEEE T. Automat. Contr. 52 (2007) 736-742.

[9] F. Plestana, Y. Shtesselb, V. Brégeaulta, A. Poznyakc, New methodologies for adaptive sliding mode control, Int. J. Control. 83 (2010) 1907-1919.

[10] Q.P. Ha, Q.H. Nguyen, D.C. Rye, H.F. Durrant-Whyte, Fuzzy Sliding-Mode Controllers with Applications, IEEE Trans. Ind. Electron. 48 (2001) 38-46.

[11] H.F. Ho, Y.K. Wong, A.B. Rad, Adaptive Fuzzy Sliding Mode Control Design: Lyapunov Approach, 5th Asian Control Conference, 2004.

[12] A. Oveisi, M. Gudarzi, Adaptive sliding mode vibration control of a nonlinear smart beam: a comparison with self-tuning Ziegler-Nichols PID controller, J. Low Freq. Noise V. A. 32 (2013) 41-62.

[13] R.J. Wai, C.M. Lin, C.F. Hsu, Adaptive fuzzy sliding-mode control for electrical servo drive, Fuzzy Set. Syst. 143 (2004) 295-310.

[14] A. Gholami, A.H.D. Markazi, A new adaptive fuzzy sliding mode observer for a class of MIMO nonlinear systems, Nonlinear Dyn. 70 (2012) 2095–2105.

[15] A. Oveisi, M. Gudarzi, M.M. Mohammadi, A. Doosthoseini, Modeling, identification and active vibration control of a funnel-shaped structure used in MRI throat, J. Vibroeng. 15 (2013) 1392-8716.

[16] S.M. Hasheminejad, A.H. Rabiee, M. Jarrahi, A.H.D. Markazi, Active vortex-induced vibration control of a circular cylinder at low Reynolds numbers using an adaptive fuzzy sliding mode controller, J. Fluid. Struc. 50 (2014) 49-65.

[17] M. Gudarzi, A. Oveisi, M.M. Mohammadi, Robust active vibration control of a rectangular piezoelectric laminate flexible thin plate: An LMI-based approach, International Review of Mechanical Engineering 6 (2012) 1217-1227.

[18] L.X. Wang, A Course in Fuzzy Systems and Control, Prentice Hall, Englewood Cliffs, 1997.

[19] S.J. Qin, An overview of subspace identification, Comput. Chem. Eng. 30 (2006) 1502-1513.

[20] S.M. Hasheminejad, A. Oveisi, Active vibration control of an arbitrary thick smart cylindrical panel with optimally placed piezoelectric sensor/actuator pairs, Int. J. Mech. Mater. Des. (2015).

[21] A. Oveisi, T. Nestorović, Robust mixed H2/H8 active vibration controller in attenuation of smart beam, F. U. Mech. Eng. 12 (2014) 235-249.

[22] T. Nestorović, M. Trajkov, Optimal actuator and sensor placement based on balanced reduced models, Mech. Syst. Signal Pr. 36 (2013) 271-289.

[23] T. Nestorović, M. Trajkov, S. Garmabi, Optimal placement of piezoelectric actuators and sensors on a smart beam and a smart plate using multi-objective genetic algorithm, Smart Struct. Syst. 15 (2015) 1041-1062.

[24] Scilab Enterprises, Scilab: Free and Open Source software for numerical computation, (2012).

[25] S. Boyd, L. E. Ghaoui, E. Feron, and V. Balakrishan, Linear matrix inequality in systems and control theory, special for industrial and applied mathematics, SIAM, Philadelphia, (1994).
