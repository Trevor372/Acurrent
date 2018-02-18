#include <Acurrent.h>
#include <math.h>

// Model Functions

static inline double
m_inf(double V)
{
    //minf(v)=1./(1.+exp(-(v+thetam)/sm))
    return 1.0/(1.0+exp(-(V+30.0)/15.0));
}

static inline double 
tau_b(double V)
{ 
 // taub0+taub1/(1+exp((v+thb)/sb))
 return 10+200/(1.0+exp((V+80/10)));
}

static inline double
n_inf(double V)
{
    //ninf(v)=1./(1.+exp((v-thetan)/sn))
    return 1.0/(1.0+exp((V + 32.0)/-8.0));
}

static inline double
tau_n(double V)
{
    //taun0+taun1/(1+exp((v+thn)/sigman))
    return 1.0+100.0/(1.0+exp((V+80.0)/26.0));
}

static inline double
a_inf(double V)
{
    //ainf(v)=1/(1+exp(-(v-thetaa)/sigmaa))
    return 1.0/(1.0+exp(-(V + 50)/20.0));
}

static inline double
b_inf(double V)
{
    //binf(v)=1/(1+exp(-(v-thetab)/sigmab))
    return 1.0/(1.0+exp(-(V + 70.0)/-6.0));
}

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
    return new Acurrent();
}

static DefaultGUIModel::variable_t vars[] =
{
    { "Vm", "Membrane Potential (V)", DefaultGUIModel::OUTPUT, },
    { "Waveform", "Total Output", DefaultGUIModel::OUTPUT, },
    { "IA_INPUT", "Input from GA_Calc (A/cm^2)", DefaultGUIModel::INPUT, },
    { "CURRENT_STEP", "Input from protocol (A/cm^2)", DefaultGUIModel::INPUT, },
    { "Iapp (uA/cm^2)", "Applied Current (uA/cm^2)",
        DefaultGUIModel::PARAMETER, },
    { "V0 (mV)", "Initial membrane potential (mV)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "phi (uF/cm^2)", "Specific membrane capacitance (uF/cm^2)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_Na_max (mS/cm^2)", "Maximum Na+ conductance density (mS/cm^2)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_Na (mV)", "Sodium reversal potential (mV)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_K_max (mS/cm^2)",
        "Maximum delayed rectifier conductance density (mS/cm^2)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_K (mV)", "K+ reversal potential (mV)", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::DOUBLE, },
    { "G_A_max (mS/cm^2)",
        "Maximum transient A-type K+ conductance density (mS/cm^2)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_L (mS/cm^2)", "Maximum leak conductance density mS/cm^2",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_L (mV)", "Leak reversal potential (mV)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Rate (Hz)", "Rate of integration (Hz)", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::UINTEGER, },
    { "Va", "Variable 1 a that effects how a is calculated", 
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Vb", "Variable b that effects how b is calculated",
	DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "n", "Potassium Activation", DefaultGUIModel::STATE, },
    { "a", "A-type Potassium Activation", DefaultGUIModel::STATE, },
    { "b", "A-type Potassium Inactivation", DefaultGUIModel::STATE, },
    { "IA", "A-type Potassium Current", DefaultGUIModel::STATE, },
    { "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, };

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

/*
* Macros for making the code below a little bit cleaner.
*/

#define V (y[0])
#define n (y[1])
#define b (y[2])
#define dV (dydt[0])
#define dn (dydt[1])
#define db (dydt[2])
#define G_K  (G_K_max*n*n*n*n)
#define G_A  (G_A_max*a*a*a*b)

Acurrent::Acurrent(void) : DefaultGUIModel("Acurrent", ::vars, ::num_vars) {
    setWhatsThis(
    "<p><b>Acurrent:</b><br>This module simulates a Acurrent model neuron.</p>");
    createGUI(vars, num_vars);
    initParameters();
    update( INIT );
    refresh();
    resizeMe();
}

Acurrent::~Acurrent(void) {}

void Acurrent::execute(void) {
    systime = count * period; // time in seconds
    for (int i = 0; i < steps; ++i)
        solve(period / steps, y); // period in s
    output(0) = V * 1e-3; // convert to V
    output(1) = 1e-9*(input(0)+input(1));
    count++;
}

void Acurrent::update(DefaultGUIModel::update_flags_t flag) {
    switch (flag) {
        case INIT:
            setParameter("V0 (mV)", QString::number(V0)); 
            setParameter("phi (uF/cm^2)", QString::number(phi)); 
            setParameter("G_Na_max (mS/cm^2)", QString::number(G_Na_max)); 
            setParameter("E_Na (mV)", QString::number(E_Na)); 
            setParameter("G_K_max (mS/cm^2)", QString::number(G_K_max)); 
            setParameter("E_K (mV)", QString::number(E_K)); 
            setParameter("G_A_max (mS/cm^2)", QString::number(G_A_max)); 
            setParameter("G_L (mS/cm^2)", QString::number(G_L)); 
            setParameter("E_L (mV)", QString::number(E_L)); 
            setParameter("Iapp (uA/cm^2)", QString::number(Iapp)); 
            setParameter("Rate (Hz)", rate);
            setParameter("Va", Va);
	    setParameter("Vb", Vb);
            setState("n", n);
            setState("a", a);
            setState("b", b);
            setState("IA",IA);
            setState("Time (s)", systime);
            break;

        case MODIFY:
            V0 = getParameter("V0 (mV)").toDouble();
            phi = getParameter("phi (uF/cm^2)").toDouble();
            G_Na_max = getParameter("G_Na_max (mS/cm^2)").toDouble();
            E_Na = getParameter("E_Na (mV)").toDouble();
            G_K_max = getParameter("G_K_max (mS/cm^2)").toDouble();
            E_K = getParameter("E_K (mV)").toDouble();
            G_A_max = getParameter("G_A_max (mS/cm^2)").toDouble();
            G_L = getParameter("G_L (mS/cm^2)").toDouble();
            E_L = getParameter("E_L (mV)").toDouble();
            Iapp = getParameter("Iapp (uA/cm^2)").toDouble(); // displayed in uA/cm^2, calculated in uA/mm^2
            rate = getParameter("Rate (Hz)").toDouble();
            Va = getParameter("Va").toDouble();
	    Vb = getParameter("Vb").toDouble();
            steps = static_cast<int> (ceil(period * rate));
            /*V = V0;
            n = n_inf(V0);
            a = a_inf(V0);
            b = b_inf(V0);*/
            break;

        case PERIOD:
            period = RT::System::getInstance()->getPeriod() * 1e-6; // time in seconds
            steps = static_cast<int> (ceil(period * rate));
            break;

        default:
            break;
    }
}

void Acurrent::initParameters() {
    V0 = -60.0; // mV
    G_Na_max = 37.0; 
    G_K_max = 45.0;
    G_L = 1.0;
    G_A_max = 40;
    E_Na = 55.0; // mV
    E_K = -80;
    E_L = -70.0;
    phi = .75;
    Iapp = 20; // 1 Hz spiking
    rate = 400;
    V = V0;
    n = 0;
    Va = 50;
    Vb = 70;
    a = a_inf(V0);
    b = .47;
    count = 0;
    systime = 0;
    period = RT::System::getInstance()->getPeriod() * 1e-6; // s
    steps = static_cast<int> (ceil(period * rate)); // calculate how many integrations to perform per execution step
}

void Acurrent::solve(double dt, double *y) {
    double dydt[3];
    derivs(y, dydt);
    for (size_t i = 0; i < 3; ++i){
        y[i] += dt * dydt[i];
    }
}

void Acurrent::derivs(double *y, double *dydt) {

    a = 1.0/(1.0+exp(-(V + Va)/20.0));
    INA = G_Na_max*pow(m_inf(V),3)*(1-n)*(V-E_Na);
    IA = G_A * (V - E_K); // A


    dV = -(G_L*(V-E_L) + INA + G_K*(V-E_K) + IA)+Iapp+input(0)+input(1);
    dn = phi*(n_inf(V) - n) / tau_n(V);
    db = ((1.0/(1.0+exp(-(V + Vb)/-6.0))) - b) / tau_b(V);
}
