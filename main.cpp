#include <bits/stdc++.h>

using namespace std;

//LINEAR EQUATIONS

vector<vector<double>> gaussian(vector<vector<double>> mat)
{
    int n = mat.size();
    vector<double> piv;
    for (int i=0; i<n; i++)
    {
        for (int k=i+1; k<n; k++)
        {
            if (abs(mat[i][i]) < abs(mat[k][i]))
            {
                swap(mat[i], mat[k]);
            }
        }
        double pivot=mat[i][i];
        if (pivot!=0)
        {
            piv.push_back(pivot);
            for (int j=i; j<=n; j++)
            {
                mat[i][j] /= pivot;
            }
        }
        for (int k=i+1; k<n; k++)
        {
            double factor = mat[k][i];
            for (int j=i; j<=n; j++)
            {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }
    // for(int i=0;i<n;i++)
    //     {
    //         for(int j=i;j<n+1;j++)
    //         {
    //             mat[i][j]*=piv[i];
    //         }
    //     }
    return mat;
}
vector<vector<double>> gauss_jordan(vector<vector<double>> mat)
{
    int n = mat.size();
    for (int i=0; i<n; i++)
        {
        for (int k=i+1; k<n; k++)
        {
            if (abs(mat[i][i]) < abs(mat[k][i]))
            {
                swap(mat[i], mat[k]);
            }
        }
        double pivot=mat[i][i];
        if (pivot!=0)
        {
            for (int j=i; j<=n; j++)
            {
                mat[i][j]/=pivot;
            }
        }
        for (int k=i+1; k<n; k++)
        {
            double factor=mat[k][i];
            for (int j=i; j<=n; j++)
            {
                mat[k][j]-=factor*mat[i][j];
            }
        }
    }
    for (int i=n-1; i>=0; i--)
    {
        for (int k=i-1; k>=0; k--)
        {
            double factor=mat[k][i];
            for (int j=i;j<=n; j++)
            {
                mat[k][j]-=factor*mat[i][j];
            }
        }
    }
    return mat;
}

//NON-LINEAR EQUATIONS
double f(vector<double> vf, double x)
{
    double a=vf[0];
    double b=vf[1];
    double c=vf[2];
    double d=vf[3];
    double e=vf[4];
    return (a*x*x*x*x)+(b*x*x*x)+(c*x*x)+(d*x)+e;
}
double fp(vector<double> vf, double x)
{
    double a=vf[0];
    double b=vf[1];
    double c=vf[2];
    double d=vf[3];
    double e=vf[4];
    return ((4*a*x*x*x)+(3*b*x*x)+(2*c*x)+d);
}

vector<double> newtonraphson(vector<double> vf,double x0)
{
    vector<double> vc;
    vc.push_back(x0);
    double fx=f(vf,x0);
    double fpx=fp(vf,x0);
    if(fpx==0){return vc;}
    else{
        double x1=x0-(fx/fpx);
        vc.push_back(x1);
        while(abs(x1-x0)>0.0001)
        {
            x0=x1;
            fx=f(vf,x0);
            fpx=fp(vf,x0);
            if(fpx==0){return vc;}
            else{
                x1=x0-(fx/fpx);
                vc.push_back(x1);
            }
        }
        return vc;
    }
}
vector<double> secant(vector<double> vf,double x1, double x2)
{
    double fx1=f(vf,x1);
    double fx2=f(vf, x2);
    vector<double> vc;
    double x3=((x1*fx2) - (x2*fx1))/(fx2-fx1);
    double fx3=f(vf, x3);
    vc.push_back(x3);
    while(abs(fx3)>0.0001)
    {
        x1=x2;
        x2=x3;
        fx1=f(vf, x1);
        fx2=f(vf, x2);

        if(fx1==fx2){return vc;}
        x3=((x1*fx2) - (x2*fx1))/(fx2-fx1);
        vc.push_back(x3);
        fx3=f(vf, x3);
    }
    return vc;
}

//DIFFERENTIAL EQUATIONS SOLVE
double f(double x,double y,double a, double b,  int type)
{
    if(type==1)
    {
        return a*sin(x);
    }
    if(type==2)
    {
        return a*x + b*y;
    }
    
}

pair<double,double> rungekutta(double h,double x,double y, double a, double b, int t)
{
    double k1= (h * f(x,y,a,b,t));
    double k2= (h* f((x + (h/2)),(y + (k1/2)),a,b,t));
    double k3= (h* f((x + (h/2)),(y + (k2/2)),a,b,t));
    double k4= (h* f((x + h),(y + k3),a,b,t));

    y= y + ((k1 + (2*k2) + (2*k3) + k4)/6.0);
    x= x + h;
    pair<double,double> xy= {x,y};
    return xy;
}

void show_rk(double a, double b, int t)
{
    vector<pair<double,double> > vxy;
    pair<double,double> xy={0.0,0.0};
    for(double h=0.1;xy.first<=(4*3.14);)
    {
        cout<<"x = "<<xy.first<<", y = "<<xy.second<<endl;
        vxy.push_back(xy);
        xy=rungekutta(h, xy.first, xy.second,a,b,t);
    }
}
//MATRIX INVERSION


//show_matrix
void show_matrix(vector<vector<double>> mat)
{
    int r=mat.size();
    int c=mat[0].size();
    cout<<"RESULTING MATRIX:"<<endl;
    for(int i=0;i<r;i++)
    {
        for(int j=0;j<c;j++)
        {
            cout<<mat[i][j]<<" ";
        }cout<<endl;
    }
}

//prompt

void returnToMainMenu() {
    cout << "\nPress Enter to return to the main menu...";
    cin.ignore();
    cin.get();
    system("CLS");
}
void showMainMenu() {
    cout<<"NUMERICAL OPERATIONS: "<<endl;
    cout<<"1. LINEAR EQUATIONS"<<endl;
    cout<<"2. NON-LINEAR EQUATIONS"<<endl;
    cout<<"3. DIFFERENTIAL EQUATIONS"<<endl;
    cout<<"4. MATRIX INVERSION"<<endl;
    cout<<"0. exit"<<endl;
    cout<<"Enter your choice: (1/2/3/4/0) ";
}
void prompt()
{
    showMainMenu();
    int x;
    cin>>x;
    if(x==1)
    {
        cout<<"Enter number of variables: ";
        int v;
        cin>>v;
        vector<vector<double>> mat(v,vector<double>(v+1));
        cout<<"Enter the coefficients: "<<endl;
        for(int i=0;i<v;i++)
        {
            for(int j=0;j<v+1;j++)
            {
                double cf;
                cin>>cf;
                mat[i][j]=cf;
            }
        }
        cout<<"Choose Numerical Methods: "<<endl;
        cout<<"1. Jacobi Iterative"<<endl;
        cout<<"2. Gauss-Seidel Iterative"<<endl;
        cout<<"3. Gauss Elimination"<<endl;
        cout<<"4. Gauss-Jordan Elimination"<<endl;
        cout<<"5. LU factorization"<<endl;
        cout<<"Enter your choice: (1/2/3/4/5) ";
        int c;
        cin>>c;
        if(c==1){/*jacobi*/}
        if(c==2){/*gauss-seidel*/}
        if(c==3)//gaussian
        {
            show_matrix(gaussian(mat));
            vector<vector<double>> gm=gauss_jordan(mat);
            cout<<"Result:"<<endl;
            for(int i=0;i<gm.size();i++)
            {
                cout<<"x"<<i+1<<"="<<gm[i][gm.size()]<<endl;
            }
        }
        if(c==4)//gauss-jordan
        {
            show_matrix(gauss_jordan(mat));
            vector<vector<double>> gm=gauss_jordan(mat);
            cout<<"Result:"<<endl;
            for(int i=0;i<gm.size();i++)
            {
                cout<<"x"<<i+1<<"="<<gm[i][gm.size()]<<endl;
            }
        }
        returnToMainMenu();
    }
    if(x==2)
    {
        cout<<"Enter degree of equation: (1-4) ";
        int deg;
        cin>>deg;
        deg++;
        vector<double> vec(5,0);
        cout<<"Enter the coefficients: ";
        for(int i=(5-deg);i<5;i++)// Took as many variables as we need, others will be 0 in the vector 'vec'
        {
            double cf;
            cin>>cf;
            vec[i]=cf;
        }
        cout<<"Choose Numerical Methods: "<<endl;
        cout<<"1. Bisection"<<endl;
        cout<<"2. False-Position"<<endl;
        cout<<"3. Secant"<<endl;
        cout<<"4. Newton-Raphson"<<endl;
        cout<<"Enter your choice: (1/2/3/4) ";
        int choice;
        cin>>choice;
        if(choice==1){/*add bisection code here*/}
        if(choice==2){/*add false-position code here*/}
        if(choice==3)//secant
        {
            vector<double> vc=secant(vec,0,1);
            cout<<"x = "<<vc[vc.size()-1]<<endl;
            cout<<"Iterations taken: "<<vc.size()<<endl;
        }
        if(choice==4)//newton-raphson
        {
            vector<double> vc=newtonraphson(vec,0);
            cout<<"x = "<<vc[vc.size()-1]<<endl;
            cout<<"Iterations taken: "<<vc.size()<<endl;
        }
        returnToMainMenu();
    }
    if(x==3)
    {
        cout<<"Type of differential equation to solve: "<<endl;
        cout<<"1. dy/dx = asin(x)"<<endl;
        cout<<"2. dy/dx = ax + by"<<endl;
        int t;cin>>t;
        double a,b;
        if(t==1){cout<<"a = ";cin>>a;b=0;}
        if(t==2){cout<<"a = ";cin>>a;cout<<"b = ";cin>>b;}

        show_rk(a,b,t);
        returnToMainMenu();
    }
    if(x==4)
    {
        //matrix inversion
        returnToMainMenu();
    }
    if(x==0)
    {
        cout << "Exiting program." << endl;
        returnToMainMenu();
    }
}

int main()
{
    while(1)
    {
        prompt();
    }
}
