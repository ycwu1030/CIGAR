// Try to follow VEGAS's algorithm to integrate over (1+x^2)/(1-x^2) from very close to -1 to very close to 1

double Integrand(double x)
{
    return (1.0+x*x)/(1.0-x*x);
}

double VEGAS1D(double xmin, double xmax, double &error)
{
    // Constants declear
    const int N = 80 + 1;  // The fixed number of increments from xmin to xmax
    const int K = 1000; // The number used to refined the increments


    // Initiation
    double xi[N];
    for (int i = 0; i < N; i++)
    {
        xi[i] = xmin + (xmax-xmin)/(N - 1.0)*i;
    }
    return 0;
    

}
