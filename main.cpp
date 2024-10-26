#include <bits/stdc++.h>

using namespace std;

//LINEAR EQUATIONS

//Checking diagonally dominant for jacobi and gauss seidel
bool diagonally_dominant(vector<vector<double>>&A) 
{
    int n=A.size();
    for (int i=0;i<n;i++) 
    {
        double sum=0.0;
        for(int j=0;j<n;j++) 
        {
            if (i!=j) 
            {
                sum+=(A[i][j]);
            }
        }
        if (fabs(A[i][i])<sum) 
        {
            return false;
        }
    }
    return true;
}
//Jacobi iteration function
vector<double>Jacobi_iterative(vector<vector<double>>& mat,vector<double>&X,int Max_Ite,double tol) 
{
    int n=mat.size();
    vector<vector<double>>A(n, vector<double>(n));
    vector<double>B(n);
    for (int i=0; i<n;i++) 
    {
        for (int j=0;j<n;j++) 
        {
            A[i][j]=mat[i][j];
        }
        B[i]=mat[i][n];
    }
    vector<double>res(n, 0.0);

    if (!diagonally_dominant(A)) 
    {
        cout << "Matrix isn't diagonally dominant. Jacobi Iteration is not possible." << endl;
        cout<<"May produce wrong results "<<endl<<endl;
    }

    vector<double>preX=X;
    for (int iter=1; iter<=Max_Ite;iter++) 
    {
        for (int i=0;i<n;i++) 
        {
            double sum=B[i];
            for (int j=0;j<n;j++)
             {
                if (j!=i) 
                {
                    sum-=A[i][j]*preX[j];
                }
            }
            X[i]=sum/A[i][i];
        }

        double diff=0.0;
        for (int i=0;i<n;i++)
         {
            diff=max(diff,fabs(X[i]-preX[i]));
        }

        if (diff<tol)
         {cout<<endl;
            cout<<"Converged in "<<iter<< " iterations."<<endl;
            return X;
        }
        preX=X;
    }
    cout<<endl;
        cout <<"Did not converge in " <<Max_Ite <<" iterations." <<endl;
        cout<<endl;
    return X;
}

//GaussSeidel iteration function 
vector<double> GaussSeidel_iterative(vector<vector<double>>& mat, vector<double>& X, int Max_Ite, double tol)
{
    int n=mat.size();
    vector<vector<double>>A(n, vector<double>(n));
    vector<double>B(n);
    
    for (int i=0;i<n;i++)
     {
        for (int j=0;j<n;j++) 
        {
            A[i][j]=mat[i][j];
        }
        B[i]=mat[i][n];
    }
    cout<<endl;

    if (!diagonally_dominant(A)) 
    {
        cout<<"Matrix isn't diagonally dominant. Gauss-Seidel Iteration may not converge."<<endl;
        cout<<endl<<"May create wrong results"<<endl<<endl;
    }

    vector<double>preX=X;
    for(int iter=1;iter<=Max_Ite;iter++)
     {
        for(int i=0;i<n;i++)
         {
            double sum=B[i];
            for (int j=0; j<n;j++)
             {
                if (j!=i) 
                {
                    sum-=A[i][j]*X[j];
                }
            }
            X[i]=sum/A[i][i];
        }

        double diff=0.0;
        for (int i=0;i<n;i++) 
        {
            diff=max(diff,fabs(X[i]-preX[i]));
        }

        if (diff<tol) 
        {
            cout<<endl;
            cout<<"Converged in "<<iter<<" iterations."<<endl;
            return X;
        }
        preX=X;
    }
    cout<<endl;
    cout <<"Did not converge in " <<Max_Ite <<" iterations." <<endl;

    return X;
}


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

bool isUniqueRoot(vector<double> vc, double root)
{
    for(double r:vc)
    {
        if(abs(r-root)<0.0001){return false;}
    }
    return true;
}
vector<double> newtonraphson(vector<double> vf,int maxAt=50)
{
    vector<double> vc;
    double tolr=0.001;

    for(int i=-maxAt;i<maxAt;i++)
    {
        double x0= i - maxAt / 2;
        double x1;
        int maxIt=100;

        for(int it=0;it<maxIt;it++)
        {
            double fx=f(vf,x0);
            double fpx=fp(vf,x0);

            if(fabs(fpx)<tolr){break;}

            x1= x0- fx/fpx;

            if(fabs(x1-x0)<tolr)
            {
                if(isUniqueRoot(vc,x1))
                {
                    vc.push_back(x1);
                }
                break;
            }
            x0=x1;
        }
    }
    sort(vc.begin(),vc.end());
    return vc;
}
vector<double> secant(vector<double> vf, int maxAt=50)
{
    vector<double> vc;
    double tol=0.0001;

    for(int i=-maxAt;i<maxAt;i++)
    {
        double x1=i;
        double x2=i+1;
        double fx1=f(vf,x1);
        double fx2=f(vf,x2);
        int maxIt=100;
        for(int it=0;it<maxIt;it++)
        {
            if(fabs(fx2-fx1)<tol){break;}
            double x3=(x1*fx2 - x2*fx1)/(fx2-fx1);
            
            if(fabs(f(vf,x3))<tol && isUniqueRoot(vc,x3))
            {
                vc.push_back(x3);
                break;
            }
            x1=x2;
            fx1=fx2;
            x2=x3;
            fx2=f(vf,x3);
        }
    }
    sort(vc.begin(),vc.end());
    return vc;
}


//bisection method


// function of checking probable intervals of roots existence
vector<pair<double, double>>intervals(vector<double>&cof,double start,double end,double step) 
{
    vector<pair<double, double>>iv;
    double a=start;
    double b=a+step;
    while (b<=end) 
    {
        if (f(cof,a)*f(cof,b)<=0) 
        {
            iv.push_back({a, b});
        }
        a=b;
        b=a+step;
    }
    return iv;
}
//bisection algorithm funtion
double bisection(vector<double>&cof,double a,double b,double tol) 
{
    double fa=f(cof,a);
    double fb =f(cof,b);
    if (fa*fb>=0) 
    {
        cout <<"Invalid interval"<<endl;
        return NAN;
    }

    double c=a;
    while ((b-a)>=tol) 
    {
        c = (a+b)/2;
        double fc=f(cof,c);

        if (fabs(fc)<tol) 
        {
            break;
        } 
        else if (fa*fc<0)
         {
            b=c;
            fb=fc;
        } 
        else 
        {
            a=c;
            fa=fc;
        }
    }
    return c;
}
//function for finding all roots of polynomial by bisection
vector<double> bisection_method_roots( vector<double>&cof,double a,double b,double tol,double step) 
{
    vector<double>roots;
    auto iv=intervals(cof,a,b,step);

    for (const auto&intv:iv) 
    {
        double a=intv.first;
        double b=intv.second;
        double root=bisection(cof,a,b,tol);
        if (!isnan(root)) 
        {
            roots.push_back(root);
        }
    }
    return roots;}


//False Position method
vector<double>false_position_method( vector<double>& cof,double a,double b,double tol,double step) 
{vector<double>roots;
     vector<pair<double, double>>iv=intervals(cof,a,b,step);
   for(pair<double, double>intv:iv)
 {
    double a=intv.first,b=intv.second,c;
        double f_a=f(cof,a);
        double f_b=f(cof,b);
        if(f_a*f_b>=0)
            continue;
        while((b-a)>=tol)
        {
            double root=b-f_b*(b-a)/(f_b-f_a);
            double f_root=f(cof,root);
            if(fabs(f_root)<tol)
            {
                roots.push_back(root);
                break;
            }
            else if(f_a*f_root<0)
            {
                b=root;
                f_b=f_root;
            }
            else
            {
                a=root;
                f_a=f_root;
            }
        }
    }
    return roots;
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
        if(c==1)//jacobi
        {
            cout<<"Enter number of Iterations you want : "<<endl;
            int it;
            cin>>it;
            vector<double>X(v,0.0);
            vector<double>solve=Jacobi_iterative(mat,X,it,0.0001);
            cout<<"Result:"<<endl;
            for(int i=0;i<solve.size();i++)
            {
                cout<<"x["<<i+1<<"] = "<<solve[i]<<endl;
            }
            cout<<endl;
        }
        if(c==2)//gauss-seidel
        {
        /*gauss-seidel*/
            cout<<"Enter number of Iterations you want : "<<endl;
                    int it;
                    cin>>it;
            vector<double>X(v,0.0);
            vector<double>solve=GaussSeidel_iterative(mat,X,it,0.0001);
              cout<<"Result:"<<endl;
            for(int i=0;i<solve.size();i++)
            {
                cout<<"x["<<i+1<<"] = "<<solve[i]<<endl;
            }
            cout<<endl;
        }
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

        if(choice==1)//bisection
        {   
            
            vector<double>vv=bisection_method_roots(vec,-50,50,0.0001,0.1);
            cout<<endl;
            cout<<"x = ";
            for(double  i : vv)
            {
                cout<<i<<" ";
            }
            cout<<endl;
        }
        if(choice==2)//false-position
        {
            vector<double>vvv=false_position_method(vec,-50,50,0.0001,0.1);
            cout<<endl;
            cout<<"x = ";
            for(double  i : vvv)
            {
                cout<<i<<" ";
            }
            cout<<endl;
        }
        if(choice==3)//secant
        {
            vector<double> vc=secant(vec);
            cout<<endl;
            cout<<"x = ";
            for(int i=0;i<vc.size();i++)
            {
                if(i==vc.size()-1)
                {
                    cout<<vc[i];
                    continue;
                }
                cout<<vc[i]<<", ";
            }cout<<endl;
        }
        if(choice==4)//newton-raphson
        {
            vector<double> vc=newtonraphson(vec);
            cout<<endl;
            cout<<"x = ";
            for(int i=0;i<vc.size();i++)
            {
                if(i==vc.size()-1)
                {
                    cout<<vc[i];
                    continue;
                }
                cout<<vc[i]<<", ";
            }cout<<endl;
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
