#include <default_gui_model.h>

class Acurrent : public DefaultGUIModel
{
	
	public:
	
		Acurrent(void);
		virtual~Acurrent(void);
	
		void execute(void);
	
	protected:
	
		void update(DefaultGUIModel::update_flags_t);
	
	private:
	
		void derivs(double *, double *);
		void solve(double, double *);
		void initParameters();
	
		double y[3];
		double period;
		int steps;
	
		double V0;
		double G_Na_max;
		double E_Na;
		double G_K_max;
		double E_K;
		double G_L;
		double E_L;
		double G_A_max;
		double Iapp;
		double a;
                double Va;
                double Vb;
		double IA;
		double INA;
		double phi;
		double rate;
		double systime;
		long long count;
};
