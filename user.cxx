//////////////////////////////////////////////////////////////////////////
////////////////        AYS model.cxx         /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:        Optimal control of the AYS model  ////////////////
//////// Last modified: 19 September 2020                 ////////////////
//////////////////////////////////////////////////////////////////////////


#include "psopt.h"

double Ac = 240; 
double Yc = 7 * pow(10 , 13); 
double Sc = 5 * pow(10 ,11); 
double phi = 4.7 * pow(10 , 10) ;
double theta = 8.57 * pow( 10 , -5) ;
double epsilon = 147 ;
double TA = 50 ;
double TS = 50 ;
double rho = 2 ;
double sigma0 = 4 * pow(10, 12); 
double beta0 = 0.03 ; 
double aupper = 0.5897;//(345 / (345 + 240));
double ylower = 0.3636;//( 4 / 11 );
//////////////////////////////////////////////////////////////////////////
/////////////////// Define auxiliary function ////////////////////////
//////////////////////////////////////////////////////////////////////////
adouble GAMMA( adouble* states, adouble* controls)
{
    adouble sigma = controls[ CINDEX(1) ];
    adouble s = states[ CINDEX(3) ];

    return pow((1 - s),rho) / ( pow((1-s) , rho) + pow((Sc * s / sigma),rho));
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{   
    adouble af = final_states[ CINDEX(1) ];
    adouble yf = final_states[ CINDEX(2) ];
    adouble L;
    L = - ( pow((af - aupper),2) + pow((yf - ylower),2) );
    // return L;//to maximize final distance from boundaries
    // return (L + tf);//to maximize final distance from boundaries & controlling time
    // return tf; // to minimize controlling time

    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    adouble sigma = controls[ CINDEX(1) ];
    adouble beta = controls[ CINDEX(2) ];

    adouble a = states[ CINDEX(1) ];
    adouble y = states[ CINDEX(2) ];

    adouble C_control;
    adouble C_boundaries;
    C_control = pow((sigma/sigma0 - 1) , 2) + pow((beta/beta0 - 1) , 2);
    C_boundaries = - ( pow((a - aupper),2) + pow((y - ylower),2) ); 
    return  C_control; // To minimize the total amount of control
    // return  C_boundaries;// To maximize distances from the boundaries in each step
    // return  (C_control + C_boundaries);// To maximize distances from the boundaries & to minimize the total amount of control
    // return 0.0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble adot, ydot, sdot;

   adouble a = states[ CINDEX(1) ];
   adouble y = states[ CINDEX(2) ];
   adouble s = states[ CINDEX(3) ];

   adouble sigma = controls[ CINDEX(1) ];
   adouble beta = controls[ CINDEX(2) ];

   adot = Yc * GAMMA(states , controls) * pow((1-a),2) * y /(phi * epsilon * Ac * (1-y)) - a*(1-a)/TA;
   ydot = y * (1-y) * (beta * theta * Ac * a) / (1-a);
   sdot = (1 - GAMMA(states , controls)) * Yc * pow((1-s),2) * y / ( epsilon * Sc * (1-y)) - s * (1-s) / TS ;

   derivatives[ CINDEX(1) ] = adot;
   derivatives[ CINDEX(2) ] = ydot;
   derivatives[ CINDEX(3) ] = sdot;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble ai = initial_states[ CINDEX(1) ];
   adouble yi = initial_states[ CINDEX(2) ];
   adouble si = initial_states[ CINDEX(3) ];

   adouble af = final_states[ CINDEX(1) ];
   adouble yf = final_states[ CINDEX(2) ];
   adouble sf = final_states[ CINDEX(3) ];


   e[ CINDEX(1) ] = ai;
   e[ CINDEX(2) ] = yi;
   e[ CINDEX(3) ] = si;

   e[ CINDEX(4) ] = af;
   e[ CINDEX(5) ] = yf;
   e[ CINDEX(6) ] = sf;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Optimal control of AYS model";
    problem.outfilename                 = "AYS.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages           = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 3;// a , y , s
    problem.phases(1).ncontrols 		= 2;// sigma , beta
    problem.phases(1).nevents   		= 6;// ai , yi , si , af , yf, sf
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes             = "[50]"; 

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double aL = 0.0;
    double aU = aupper;
    double yL = ylower;
    double yU = 1.0;
    double sL = 0.0;
    double sU = 1.0;



    double u1L = pow(10, 12);//sigma0 / 4;
    double u2L = 0.015;//beta0 / 2 ;
    double u1U = sigma0 ;
    double u2U = beta0 ;


    double ai = 0.5;
    double yi = 0.5;
    double si = 0.5;
///////////////////////shelter's bounds///////////////

    double af_lower = 0.001;
    double yf_lower = 0.3671;
    double sf_lower = 0.6834;

    double af_upper = 0.4273;
    double yf_upper = 0.8134;
    double sf_upper = 1.0;
///////////////////////////////////////
    double TGuess = 100.0;

    problem.phases(1).bounds.lower.states(1) = aL;
    problem.phases(1).bounds.lower.states(2) = yL;
    problem.phases(1).bounds.lower.states(3) = sL;


    problem.phases(1).bounds.upper.states(1) = aU;
    problem.phases(1).bounds.upper.states(2) = yU;
    problem.phases(1).bounds.upper.states(3) = sU;


    problem.phases(1).bounds.lower.controls(1) = u1L;
    problem.phases(1).bounds.lower.controls(2) = u2L;
    problem.phases(1).bounds.upper.controls(1) = u1U;
    problem.phases(1).bounds.upper.controls(2) = u2U;

    problem.phases(1).bounds.lower.events(1) = ai;
    problem.phases(1).bounds.lower.events(2) = yi;
    problem.phases(1).bounds.lower.events(3) = si;
    problem.phases(1).bounds.lower.events(4) = af_lower;
    problem.phases(1).bounds.lower.events(5) = yf_lower;
    problem.phases(1).bounds.lower.events(6) = sf_lower;

    problem.phases(1).bounds.upper.events(1) = ai;
    problem.phases(1).bounds.upper.events(2) = yi;
    problem.phases(1).bounds.upper.events(3) = si;
    problem.phases(1).bounds.upper.events(4) = af_upper;
    problem.phases(1).bounds.upper.events(5) = yf_upper;
    problem.phases(1).bounds.upper.events(6) = sf_upper;


    problem.phases(1).bounds.upper.path(1) = 0.0;
    problem.phases(1).bounds.lower.path(1) = 0.0;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 5.0;///???
    problem.phases(1).bounds.upper.EndTime      = 50.0;///?????

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae           	= &dae;
    problem.events 		    = &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			            = problem.phases(1).nodes(1);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    DMatrix x_guess    =  zeros(nstates,nnodes);

    x_guess(1,colon()) = ai*ones(1,nnodes);
    x_guess(2,colon()) = yi*ones(1,nnodes);
    x_guess(3,colon()) = si*ones(1,nnodes);

    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0,TGuess,nnodes);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.mesh_refinement             = "automatic";
    algorithm.collocation_method = "Legendre";
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");
    lambda.Save("lambda.dat");
    H.Save("H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name+": states", "time (s)", "states","a y s");

    plot(t,u(1,colon()),problem.name+": Sigma","time (s)", "controls", "Sigma");
    plot(t,u(2,colon()),problem.name+": Beta","time (s)", "controls", "Beta");


    plot(t,x,problem.name+": states", "time (s)", "states","a y s",
                             "pdf", "AYS_states.pdf");

    plot(t,u(1,colon()),problem.name+": Sigma","time (s)", "Sigma", "Sigma",
                             "pdf", "AYS_sigma.pdf");

    plot(t,u(2,colon()),problem.name+": Beta","time (s)", "Beta", "Beta",
                             "pdf", "AYS_beta.pdf");                             
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
