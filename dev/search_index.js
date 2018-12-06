var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#RationalExpectations-documentation-1",
    "page": "Home",
    "title": "RationalExpectations documentation",
    "category": "section",
    "text": "This packages provides functions that solve (linear) rational expectations model. The methods implemented at this point are gensys and Klein(2000)"
},

{
    "location": "sims/#",
    "page": "Sims",
    "title": "Sims",
    "category": "page",
    "text": ""
},

{
    "location": "sims/#RationalExpectations.sims-NTuple{4,AbstractArray}",
    "page": "Sims",
    "title": "RationalExpectations.sims",
    "category": "method",
    "text": "sims(G0,G1,Pi,Psi)\n\nThis function solves a rational expectations model using Sims(2000) (aka gensys). The model should be written in the following way:\n\nGamma_0 mathbfx_t+1 = Gamma_1 mathbfx_t + Psi mathbfvarepsilon_t + Pi mathbfta_t\n\nThe arguments of the function are:\n\nG0 is Gamma_0\nG1 is Gamma_1\nPi is Pi\nPsi is Psi\n\nWhere Gamma_0 can be singular. eta_t is the expectation error and E_t(nu_t+1)=0 by construction.\n\nSims`s solver does not require that we say what variables are jump variables.\n\nValue\n\nThe function returns theta1 and theta2, such that:\n\ny_t = Theta_1 y_t-1 + Theta_2 varepsilon_t\n\nNote\n\nAlthough the original algorithm allows cases with multiple equilibrium (sunspots), this has not been implemented thus far.\n\n\n\n\n\n"
},

{
    "location": "sims/#RationalExpectations.irf-Tuple{AbstractArray,AbstractArray,Int64,Any}",
    "page": "Sims",
    "title": "RationalExpectations.irf",
    "category": "method",
    "text": "irf(model::SimsSol,t::Int,shock)\nirf(Theta1,Theta2,t,shock)\n\nGenerates the IRF from the matrices calculate using the gensys.\n\nThis function allows you to pass the matrices directly or to provide the whole model.\n\nmodel receives the full output of the function sims\nTheta1 is the theta1 matrix from sims\nTheta2 is the theta1 matrix from sims\nt is the number of periodos to be simulated\nshock is the size of the shock\n\nSee also sims(G0,G1,Pi,Psi)\n\n\n\n\n\n"
},

{
    "location": "sims/#Sims-Solver-1",
    "page": "Sims",
    "title": "Sims Solver",
    "category": "section",
    "text": "\nsims(G0::AbstractArray,G1::AbstractArray,Pi::AbstractArray,Psi::AbstractArray)\n\nirf(Theta1::AbstractArray,Theta2::AbstractArray,t::Int,shock)"
},

{
    "location": "sims/#Important-1",
    "page": "Sims",
    "title": "Important",
    "category": "section",
    "text": "When you add an autocorrelated shock to the system, the value on the lhs must be on t+1. Otherwise, the IRF will be wrongly computed."
},

{
    "location": "sims/#Note-1",
    "page": "Sims",
    "title": "Note",
    "category": "section",
    "text": "We don\'t force Julia to use the complex Schur decomposition. If you want it, just use sims(complex(G0),complex(G1),Pi,Psi). However, you will not be able to plot (yet) simply passing the model to irf"
},

{
    "location": "klein/#",
    "page": "Klein",
    "title": "Klein",
    "category": "page",
    "text": ""
},

{
    "location": "klein/#RationalExpectations.klein-Tuple{AbstractArray,AbstractArray,AbstractArray,Int64,AbstractArray,AbstractArray,AbstractArray}",
    "page": "Klein",
    "title": "RationalExpectations.klein",
    "category": "method",
    "text": "klein(A,B,C,t,k0,shock_exp,jumps)\n\nA,B,C see bellow\nt is the number of periods to be simulated\nk0 is the initial condition for the stable part of the system.\nshock_exp is the expected value of the future shocks, usually zero.\njump is a vector that tells the position of the jump variables\n\nSolves a rational expectation models using Klein(2000) method. The model has to be written in the following way:\n\nA E_t(mathbfx_t+1) = B mathbfx_t + C mathbfvarepsilon_t\n\nWhere mathbfx_t is the vector we are interested in and varepsilon_t is a vector of random shocks. For now, the jumps variables must be the last variables in the system.\n\nTo guarantee that there is a stable solution, we need to check the Blanchard Khan conditions: the number of eigenvalues bigger than 1 is equal to the number of non predeterminated variables in the problem (e.g. variables that depend on the expectation of its future value), which are also called jump variables.\n\nThis method accepts a matrix A that is singular.\n\n\n\n\n\n"
},

{
    "location": "klein/#Klein-Solver-1",
    "page": "Klein",
    "title": "Klein Solver",
    "category": "section",
    "text": "klein(A::AbstractArray,B::AbstractArray,C::AbstractArray,t::Int,k0::AbstractArray,shock_exp::AbstractArray,jumps::AbstractArray)"
},

{
    "location": "klein/#Important-1",
    "page": "Klein",
    "title": "Important",
    "category": "section",
    "text": "When you add an autocorrelated shock to the system, the value on the lhs must be on t+1. Otherwise, the IRF will be wrongly computed."
},

{
    "location": "example/#",
    "page": "Example",
    "title": "Example",
    "category": "page",
    "text": "We will implement the model from Galí(2008), chapter 3. We will need to write two separated models, as the sintax for the Sims solver is different of the sintaxe of the Klein solver. First, let\'s use the same calibration as Galí(2008):\nsigma = 1\nphi_pi = 1.5\nphi_y = 0.5/4\nbeta = 0.99\nphi = 1\nalph = 1/3\nep = 6\ntheta = 2/3\nrho_v = 0.5\n\nTheta = (1-alph)/(1-alph+alph*ep)\nlambda = (1-theta)*((1-beta*theta)/theta)*Theta\nkappa = lambda*(sigma+(phi+alph)/(1-alph))\nThe equations are:pi_t = beta E_t(pi_t+1) + kappatildey_t\ntildey_t = -frac1sigma(i_t - E_t(pi_t+1)) + E_t(tildey_t+1)\ni_t = rho + phi_pipi_t + phi_tildey_ttildey_t + v_t\nv_t = rho_v v_t-1 + varepsilon_vWe ignore the r^n_t term as Galí(2008) does."
},

{
    "location": "example/#Klein-1",
    "page": "Example",
    "title": "Klein",
    "category": "section",
    "text": "The ordering of the variables for this model is:x_t+1 = beginpmatrix\nv_t+1\ni_t\nE_t(pi_t+1)\nE_t(tildey_t+1)\nendpmatrixNotice that we are using mathbfv_t+1 on the left hand side. Using v_mathbft will generate the wrong matrices. We will write this in a way that is consistent with the Klein method (equation 1), so the matrices areA = beginpmatrix\n1  0  0  0\n0  1  0  0\n0  0  beta  0\n0  -1  1  sigma\nendpmatrix\n\nB = beginpmatrix\nrho_v  0  0  0\n1  0  phi_pi  phi_y\n0  0  1  -kappa\n0  0  0  sigma\nendpmatrix\n\nC = beginpmatrix\n1\n0\n0\n0endpmatrixAs in Gali(2008), we will set a initial monetary poliocy shock of 0.25. This works as an initial condition for the model. As usual, we set that the expected value of the shocks are zero. We will receive the impulse response function automatically.\nusing Plots, RationalExpectations\n\nA = [[1 0 0 0];[0 1 0 0]; [0 0 beta 0]; [0 -1 1 sigma]]\nB = [[rho_v 0 0 0];[1 0 phi_pi phi_y];[0 0 1 -kappa];[0 0 0 sigma]]\n\nC = [1;0; 0; 0]\n\nk0 = [0.25;0]\n\nt=12\n\nchoque = [0.; 0;0; 0]\n\nklein_sol = klein(A,B,C,t,k0,choque,[3 4])\nWe can plot the elements of klein_sol to see the irf"
},

{
    "location": "example/#Sims-1",
    "page": "Example",
    "title": "Sims",
    "category": "section",
    "text": "Sims methods requires that we write expectations error, e.g. eta_t^pi = pi_t - E_t-1(pi_t). We can work it to obtain the following matrices:Gamma_0 = beginpmatrix\n1  0  0  0\n0  1  0  0\n0  -1  sigma  1\n0  0  0  beta\nendpmatrix\n\nGamma_1 = beginpmatrix\n rho_v  0  0  0\n1  0  phi_y  phi_pi\n0  0  sigma  0\n0  0  -kappa  1\nendpmatrix\n\nPsi = beginpmatrix\n1\n0\n0\n0\nendpmatrix\n\nPi = beginpmatrix\n0  0\nphi_y  phi_pi\nsigma  0\n-kappa  1\nendpmatrixAnd the variables are ordered as:x_t+1 = beginpmatrix\nv_t+1\ni_t\nE_t(tildey_t+1)\nE_t(pi_t+1)\nendpmatrixSee the end of this article for the whole maths of this transformation. Here is it, in Julia:\nG0 = [[1 0 0 0];[0 1 0 0];[0 -1 sigma 1];[0 0 0 beta]]\nG1 = [[rho_v 0 0 0];[1 0 phi_y phi_pi];[0 0 sigma 0];[0 0 -kappa 1]]\nPsi = [1 0 0 0]\'\nPi = [[0 phi_pi 0 1];[0 phi_y sigma -kappa]]\'\n\nsol_sims = sims(G0,G1,Pi,Psi)\n\nresul = irf(sol_sims,12,0.25)\n"
},

{
    "location": "example/#Analytical-Solution-and-comparison-1",
    "page": "Example",
    "title": "Analytical Solution and comparison",
    "category": "section",
    "text": "Last, but not least, Galí(2008) gives an analytical solution for tildey_t and pi_t They are:tildey_t = -(1-betarho_v) Lambda_v v_t\npi_t = -kappaLambda_v v_tAnd Lambda_v = frac1(1-betarho_v)sigma(1-rho_v)+phi_y+kappa(phi_pi-rho_v)Lets put the analytical solutions in Julia:\nLambda_v = 1/((1-beta*rho_v)*(sigma*(1-rho_v)+phi_y)+kappa*(phi_pi - rho_v))\ntrue_y(v) = -(1-beta*rho_v)*Lambda_v*v\ntrue_pi(v) = -kappa*Lambda_v*v\n\ntrue_path = zeros(13,3)\ninitial_shock = 0.25\nshock = zeros(13)\nshock[2] = initial_shock\n\nfor j = 2:13\n    true_path[j,1] = rho_v*true_path[j-1] + shock[j]\n    true_path[j,2] = true_y(true_path[j,1])\n    true_path[j,3] = true_pi(true_path[j,1])\nend\nWe are able to compare the analytical solution with the estimated solutions. First, the shock on v_t\nplot(true_path[2:13,1], lab = \"Analytical Solution\")\nplot!(resul[:,1], lab = \"Gensys Answer\", linestyle = :dash))\nplot!(klein_sol[:,1], lab = \"Klein Answer\", linestyle = :dot))Here is the shock in the output gap:\nplot(true_path[2:13,2], lab = \"Analytical Solution\")\nplot!(resul[:,3], lab = \"Gensys Answer\", linestyle = :dash)\nplot!(klein_sol[:,4], lab = \"Klein Answer\", linestyle = :dot)\nAnd the shock in the inflation:\nplot(4*true_path[2:13,3], lab = \"Analytical Solution\")\nplot!(4*resul[:,4], lab = \"Gensys Answer\", linestyle = :dash)\nplot!(4*klein_sol[:,3], lab = \"Klein Answer\", linestyle = :dot)\n"
},

{
    "location": "biblio/#",
    "page": "Bibliography",
    "title": "Bibliography",
    "category": "page",
    "text": "Miao, Jianjun (2014) Economic Dynamics in discrete Time, MIT Press Book\nSims, Christopher(2000) Solving Linear Rational Expectations Models You can read it here.\nKlein, Paul (2000) Using generalized Schur form to solve multivariate linear rational expectations model Journal of Economic Dynamics and Control 24(1):1405-23\nGalí, Jordi (2008) Monetary Policy, Inflation and Business Cycle: An Introduction to the New Keynesian Framework, Princeton University Press"
},

]}
