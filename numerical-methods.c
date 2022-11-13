#include <stdio.h>
#include <math.h>
#define STRING	100
#define SIZE	10
#define MORESIZE 20

void getPoly(double polynomial[2][SIZE],int	xCount){
    int		i;

    for (i = 0; i < xCount; i++){
        printf("Enter x's power: ");
        scanf("%lf",&polynomial[0][i]);

        printf("Enter x's coefficient: ");
        scanf("%lf",&polynomial[1][i]);
    }
}
// This function helps to solve polynomial equations
double func(double polynomial[2][SIZE], double x, int xCount){
    int i;
    double result = 0;
    for(i = 0; i < xCount; i++){
        result += polynomial[1][i]*pow(x,polynomial[0][i]);
    }
    return result;
}
void derivative(double polynomial[2][SIZE], double derivativedPoly[2][SIZE], int xCount){
	int	i;
	
	for(i = 0; i < xCount; i++){
		derivativedPoly[0][i] = polynomial[0][i];
		derivativedPoly[1][i] = polynomial[1][i];
	}
	
	for(i = 0; i < xCount; i++){
		if(derivativedPoly[0][i] == 0){
			derivativedPoly[1][i] = 0;
		}else{
			derivativedPoly[1][i] *= derivativedPoly[0][i];
			derivativedPoly[0][i]--;
		}
	}
}
void getMatrix(double matrix[SIZE][SIZE], int dim){
	int	i,j;
	
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			printf("Provide a number for matrix[%d][%d]: ",i+1,j+1);
			scanf("%lf",&matrix[i][j]);
		}
		printf("\n");
	}
}
void getResultMatrix(double resultMatrix[SIZE][1], int dim){
	int	i,j;
	
	for(i = 0; i < dim; i++){
		printf("Provide a number for result matrix[%d][1]: ",i+1);
		scanf("%lf",&resultMatrix[i][0]);
		printf("\n");
	}
}
void gaussMaker(double matrix[SIZE][SIZE], double resultMatrix[SIZE][1], int dim){
	int i,j,k;
	double	toZero; // Matrix element that we want to make it zero
	double	diagonal; // Number which is in the diagonal
	
	for(i = 0; i < dim; i++){
		diagonal = matrix[i][i];
		for(j = i+1; j < dim; j++){
			toZero = matrix[j][i];
			for(k = i; k < dim; k++){
				matrix[j][k] = matrix[j][k]*diagonal/toZero;
				matrix[j][k] = matrix[i][k] - matrix[j][k];
			}
			resultMatrix[j][0] = resultMatrix[j][0]*diagonal/toZero;
			resultMatrix[j][0] = resultMatrix[i][0] - resultMatrix[j][0];
		}
	}
}
void gaussSolver(double matrix[SIZE][SIZE], double resultMatrix[SIZE][1], double variables[SIZE], int dim){
	int 	i,j;
	
	for(i = dim-1; i >= 0; i--){
		variables[i] = resultMatrix[i][0]/matrix[i][i];
		for(j = dim-1; j >= 0; j--){
			resultMatrix[j][0] = resultMatrix[j][0] - matrix[j][i]*variables[i];
		}
	}
	
	for(i = 0; i < dim; i++){
		printf("x%d= %lf ",i+1 ,variables[i]);
	}
}	
double	getDeterminant(double matrix[SIZE][SIZE], int dim){
	int 	i,j,k;
	double	diagonal; // Number which is in the diagonal
	double	toZero; // Matrix element that we want to make it zero
	double	det = 1;
	double 	controlMatrix[SIZE][SIZE]; // We don't want to change our matrix, so we use stand-in matrix here
	
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			controlMatrix[i][j] = matrix[i][j];
		}
	}
	
	for(i = 0; i < dim; i++){
		diagonal = controlMatrix[i][i];
		for(j = i+1; j < dim; j++){
			toZero = controlMatrix[j][i];
			for(k = i; k < dim; k++){
				controlMatrix[j][k] = controlMatrix[j][k]*diagonal/toZero;
				controlMatrix[j][k] = controlMatrix[i][k] - controlMatrix[j][k];
			}
			det = det*toZero/diagonal; // We track changes in determinant here 
			det *= -1;				   // And then we will apply it in another for loop	
		}
	}
	
	for(i = 0; i < dim; i++){
		det *= controlMatrix[i][i];	   // Merging changes to determinant and diagonals' multiplication
	}
	
	return det;
}
void invertMatrix(double matrix[SIZE][SIZE], int dim){
	int 	i,j,k;
	double 	diagonal;
	double 	unitMatrix[SIZE][SIZE] = {{1,0,0,0,0,0,0,0,0,0},
									  {0,1,0,0,0,0,0,0,0,0},
									  {0,0,1,0,0,0,0,0,0,0},
									  {0,0,0,1,0,0,0,0,0,0},
									  {0,0,0,0,1,0,0,0,0,0},
									  {0,0,0,0,0,1,0,0,0,0},
									  {0,0,0,0,0,0,1,0,0,0},
									  {0,0,0,0,0,0,0,1,0,0},
									  {0,0,0,0,0,0,0,0,1,0},
									  {0,0,0,0,0,0,0,0,0,1}}; // We could have define it by for loops but 'ne gerek var'
							  
	
			  
	for(i = 0; i < dim; i++){
		diagonal = matrix[i][i];
		for(j = 0; j < dim; j++){
			matrix[i][j] /= diagonal; 
			unitMatrix[i][j] /= diagonal;
		}	
		for(j = 0; j < dim; j++){
			if(i != j){
				diagonal = matrix[j][i]; // We actually do not use it here as diagonal but we don't want to have a new variable for this
				for(k = 0; k < dim; k++){
					matrix[j][k] = matrix[j][k] - matrix[i][k]*diagonal;
					unitMatrix[j][k] = unitMatrix[j][k] - unitMatrix[i][k]*diagonal;
				}
			}
		}
	}
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			printf("%lf ",unitMatrix[i][j]);
		}
		printf("\n");
	}
}
double	forwardDiff(double dataValues[MORESIZE], int whichDiffer, int x){
	double zero;
	
	if(whichDiffer > 0){
		zero = forwardDiff(dataValues,whichDiffer-1,x+1) - forwardDiff(dataValues,whichDiffer-1,x);
	}else{
		zero = dataValues[x];
	}
	
	return zero;
}
double	factorial(int x){
	double i,result = 1;
	
	for(i = 1; i <= x; i++){
		result *= i;
	}
	return result;
}
int main(){
    int		choice;
    // Matrix for polynomials
    double	polynom[2][SIZE] = {0};
    double	polyDeriv[2][SIZE];
    double 	matrix[SIZE][SIZE];
    double	resultMatrix[SIZE][1];
    double	variables[SIZE] = {0}; // We'll use this in Gauss elimination method to keep variables
    double	dataValues[MORESIZE] = {0}; // We'll use this in Gregory-Newton enterpolation to keep values of datas
    int     xCount; // Number of x's in polynomial
    double	a;
    double 	b; 
    double  c = 0;
    double  d; // Starting point for Newton-Raphson
    double	x; // Numerical differentiation value and Grergory-Newton Enterpolation value
    double	h; // Used for representing gaps
    double	temp;
    double	tmp[SIZE];
	double	toVariable; // var[0] = result[] - matrix[1]*var[1] - matrix[2]var[2]... This goes on. toVariable calculates matrix*var parts
    int		iterMax;
    int 	n;
    double	tolerance;
    double	diagonal;
    double	result = 0; // Numerical differentiation and Simpson's rule result
    double	entResult = 0;
    double	initVal; // Initial value for integral
    double	finalVal; // Final value for integral
    int		i,j,k;
    double	iDbl, jDbl; // For loop elements in Simpson's rule and trapezoidal rule
    int		flag = 1;
    int		dim; // Variable for matrixes' dimension
    char	menu[11][STRING] = {"Quit: 0",
                                "Bisection: 1",
                                "Regula-Falsi: 2",
                                "Newton-Raphson: 3",
                                "Inverse Matrix: 4",
                                "Gauss Elimination: 5",
                                "Gauss-Seidel :6",
                                "Numerical Differentiation (Backward,Central,Forward Differences) : 7",
                                "Simpson's Rule: 8",
                                "Trapezoidal Rule: 9",
                                "Gregory-Newton Interpolation: 10\n"};



    printf("Hello! Please choose your calculation method.\n");
    for(i=0;i<11;++i){
        printf("%s\n",menu[i]);
    }
    printf("Choice: ");
    scanf("%d",&choice);
    printf("\n");

    if(choice == 0){
        printf("Good Bye!");
    }
    else if(choice == 1){
        printf("Enter a max iteration value: ");
        scanf("%d",&iterMax);
        printf("Enter a tolerance value: ");
        scanf("%lf",&tolerance);
        printf("Enter an a value for [a,b] (Which f(a) is smaller than 0) : ");
        scanf("%lf",&a);
        printf("Enter a b value for [a,b] (Which f(b) is bigger than 0) : ");
        scanf("%lf",&b);
        printf("Please enter how many x do you want in your polynomial (For example ax2 + bx + c has 3 x values): ");
        scanf("%d",&xCount);
        getPoly(polynom, xCount);
		
		for(i = 0; i < iterMax - 1; i++){
            c = (b + a)/2; // Defining new midpoints in each iteration
            if((b - a)/ pow(2,i) < tolerance){
               	i = iterMax;
               	flag = 0;
           	}
           	if(flag == 1){
           		if(func(polynom,c,xCount)*func(polynom,a,xCount) < 0){
               		printf("\nIteration %d: ",i+1);
               		printf("\nStart: %lf",a);
               		printf("\nMid: %lf",c);
               		printf("\nEnd: %lf",b);
               		printf("\nf(start): %lf",func(polynom,a,xCount));
               		printf("\nf(mid): %lf",func(polynom,c,xCount));
               		printf("\nf(end): %lf\n",func(polynom,b,xCount));
              		b = c;
            	}else{
            		printf("\nIteration %d: ",i+1);
               		printf("\nStart: %lf",a);
               		printf("\nMid: %lf",c);
               		printf("\nEnd: %lf",b);
               		printf("\nf(start): %lf",func(polynom,a,xCount));
               		printf("\nf(mid): %lf",func(polynom,c,xCount));
               		printf("\nf(end): %lf\n",func(polynom,b,xCount));
           			a = c;
				}
			}
        }
        printf("\nResult: %lf",c);
    }
    else if(choice == 2){
    	printf("Enter a max iteration value: ");
        scanf("%d",&iterMax);
        printf("Enter a tolerance value: ");
        scanf("%lf",&tolerance);
        printf("Enter an a value for [a,b] (Which f(a) is smaller than 0) : ");
        scanf("%lf",&a);
        printf("Enter a b value for [a,b] (Which f(b) is bigger than 0) : ");
        scanf("%lf",&b);
        printf("Please enter how many x do you want in your polynomial (For example ax2 + bx + c has 3 x values): ");
        scanf("%d",&xCount);
        getPoly(polynom, xCount);
        
        for(i = 0; i < iterMax - 1; i++){
           	c = (a*func(polynom,b,xCount) - b*func(polynom,a,xCount))/(func(polynom,b,xCount) - func(polynom,a,xCount)); // Defining new points in each iteration
           	if((b - a)/ pow(2,i) < tolerance){
               	i = iterMax;
               	flag = 0;
           	}
         	if(flag == 1){
           		if(func(polynom,c,xCount)*func(polynom,a,xCount) < 0){
               		printf("\nIteration %d: ",i+1);
               		printf("\nStart: %lf",a);
               		printf("\nEnd: %lf",b);
               		printf("\nPoint: %lf",c);
               		printf("\nf(start): %lf",func(polynom,a,xCount));
               		printf("\nf(end): %lf",func(polynom,b,xCount));
               		printf("\nf(point): %lf\n",func(polynom,c,xCount));
             		b = c;
            	}else{
            		printf("\nIteration %d: ",i+1);
               		printf("\nStart: %lf",a);
               		printf("\nEnd: %lf",b);
               		printf("\nPoint: %lf",c);
               		printf("\nf(start): %lf",func(polynom,a,xCount));
               		printf("\nf(end): %lf",func(polynom,b,xCount));
               		printf("\nf(point): %lf\n",func(polynom,c,xCount));
           			a = c;
				}
			}
       	}
       	printf("\nResult: %lf",c);
	}
	else if(choice == 3){
		printf("Enter a max iteration value: ");
        scanf("%d",&iterMax);
        printf("Enter a tolerance value: ");
        scanf("%lf",&tolerance);
        printf("Enter starting point: ");
        scanf("%lf",&d);
        printf("Please enter how many x do you want in your polynomial (For example ax2 + bx + c has 3 x values): ");
        scanf("%d",&xCount);
        
        getPoly(polynom, xCount);
        derivative(polynom,polyDeriv,xCount);
        
        for(i = 0; i < iterMax - 1; i++){
        	if(func(polyDeriv,d,xCount-1) == 0){
        		printf("Inside somewhere some values couldn't be calculated (0/0 etc.). Try to use another starting point.");
        		i = iterMax;
        		flag = 0;
			}
        	c = d - func(polynom,d,xCount)/func(polyDeriv,d,xCount-1);
            if(fabs(c - d) < tolerance){
                i = iterMax;
                flag = 0;
            }
            if(flag == 1){
            	printf("\nIteration%d: \n",i+1);
            	printf("xn: %lf",c);
            	printf("\nf(xn): %lf\n ",func(polynom,c,xCount));
            	c = c + d;
            	d = c - d;
            	c = c - d;
			}
        }
	}
	else if(choice == 4){
		printf("Please give a dimension for a square matrix: ");
		scanf("%d",&dim);
		getMatrix(matrix,dim);
		if(getDeterminant(matrix,dim) == 0){
			printf("That matrix has no invert. Because it's determinant equals 0.");
			flag = 0;
		}
		if(flag == 1){
			invertMatrix(matrix,dim);
		}
	}
	else if(choice == 5){
		printf("Please give a dimension for a square matrix (We want variable and equation counts are the same): ");
		scanf("%d",&dim);
		getMatrix(matrix,dim);
		getResultMatrix(resultMatrix,dim);
		gaussMaker(matrix,resultMatrix,dim);
		gaussSolver(matrix,resultMatrix,variables,dim);
	}
	else if(choice == 6){
		printf("Enter a max iteration value: ");
        scanf("%d",&iterMax);
        printf("Enter a tolerance value: ");
        scanf("%lf",&tolerance);
        printf("Please give a dimension for a square matrix: ");
		scanf("%d",&dim);
        printf("Please provide diagonally dominant matrix.\n");
		getMatrix(matrix,dim);
		
		for(i = 0; i < dim; i++){
			diagonal = matrix[i][i];
			for(j = 0; j < dim; j++){
				if(fabs(diagonal) < fabs(matrix[i][j])){
					flag = 0;
				}
			}
		}
		
		if(flag == 1){
			getResultMatrix(resultMatrix,dim);
			printf("Please provide starting points for variables: \n");
			for(i = 0; i < dim; i++){
				printf("x%d: ", i+1);
				scanf("%lf",&variables[i]);
			}
	
			for(i = 0; i < iterMax; i++){
				for(j = 0; j < dim; j++){
					toVariable = 0;
					for(k = 0; k < dim; k++){
						if(k != j){
							toVariable += matrix[j][k]*variables[k]; // e.g x = (6 - 7y -4z)/5 Here we take other undefined parts of equation (y and z) as constant, so we calculate it here
						}
					}
					variables[j] = (resultMatrix[j][0] - toVariable)/matrix[j][j];
				}
				for(j = 0; j < dim; j++){ // If 'if' block runs else block for once we'll understand not all variables' change are below tolerance value
					if((fabs(variables[j] - tmp[j]) < tolerance) && i > 0){
						flag = 0;
					}else{
						flag = 1;
					}
				}
				printf("Iteration%d: ",i+1);
				for(j = 0; j < dim; j++){
					printf("x%d %lf ", j+1,variables[j]);
					tmp[j] = variables[j];
				}
				printf("\n");
				if(flag == 0){
					i = iterMax;
				}
			}
		}else{
			printf("This matrix is not diagonally dominant.");
		}
	}
	else if(choice == 7){
		printf("Choose which differantiation you want: ");
		printf("\n1- Backward Difference");
		printf("\n2- Centered Difference");
		printf("\n3- Forward Difference\n");
		printf("\nChoice: ");
		scanf("%d",&choice);
		printf("\nPlease enter how many x do you want in your polynomial (For example ax2 + bx + c has 3 x values): ");
        scanf("%d",&xCount);
		getPoly(polynom,xCount);
		printf("Which x do you want to differentiate: ");
		scanf("%lf",&x);
		printf("Enter the h value: ");
		scanf("%lf",&h);
		if(choice == 1){
			result = (func(polynom,x,xCount) - func(polynom,x-h,xCount))/h;
			printf("Result: %lf",result);
		}
		else if(choice == 2){
			result = (func(polynom,x+h,xCount) - func(polynom,x-h,xCount))/(2*h);
			printf("Result: %lf",result);
		}
		else if(choice == 3){
			result = (func(polynom,x+h,xCount) - func(polynom,x,xCount))/h;
			printf("Result: %lf",result);
		}
		else{
			printf("That's not a valid choice.");
		}
				
	}
	else if(choice == 8){
		printf("Choose the formula you want: ");
		printf("\n1- Simpson's 1/3 rule");
		printf("\n2- Simpson's 3/8 rule\n");
		printf("\nChoice: ");
		scanf("%d",&choice);
		printf("\nProvide an initial value for integral: ");
		scanf("%lf",&initVal);
		printf("Provide a final value for integral: ");
		scanf("%lf",&finalVal);
		printf("Provide n value: ");
		scanf("%d",&n);
		printf("Please enter how many x do you want in your polynomial (For example ax2 + bx + c has 3 x values): ");
        scanf("%d",&xCount);
        getPoly(polynom,xCount);
        
        a = (finalVal-initVal)/n;
        
		if(choice == 1){
			for(iDbl = initVal; iDbl < finalVal; iDbl +=a){
				jDbl = iDbl+a; // iDbl == gap's initVal jDbl == gap's finalVal
				result += (jDbl-iDbl)*(func(polynom,iDbl,xCount) + 4*func(polynom,(iDbl+jDbl)/2,xCount) + func(polynom,jDbl,xCount))/6; // Simpson's 1/3 formula
			}
			printf("%lf",result);
		}
		else if(choice == 2){
			for(iDbl = initVal; iDbl < finalVal; iDbl += a){
				jDbl = iDbl+a; // iDbl == gap's initVal jDbl == gap's finalVal
				result += (jDbl-iDbl)*(func(polynom,iDbl,xCount) + 3*func(polynom,iDbl + (jDbl-iDbl)/3,xCount) + 3*func(polynom,iDbl + 2*(jDbl-iDbl)/3,xCount) + func(polynom,jDbl,xCount))/8; // Simpson's 3/8 formula
			}
			printf("%lf",result);
		}else{
			printf("That's not a valid choice.");
		}
	}
	else if(choice == 9){
		printf("Provide an initial value for integral: ");
		scanf("%lf",&initVal);
		printf("Provide a final value for integral: ");
		scanf("%lf",&finalVal);
		printf("Provide n value: ");
		scanf("%d",&n);
		printf("Please enter how many x do you want in your polynomial (For example ax2 + bx + c has 3 x values): ");
        scanf("%d",&xCount);
        getPoly(polynom,xCount);
        
         a = (finalVal-initVal)/n;
         
         for(iDbl = initVal; iDbl < finalVal; iDbl += a){
         	jDbl = iDbl+a;
         	result += (jDbl-iDbl)*(func(polynom,iDbl,xCount) + func(polynom,jDbl,xCount))/2; // Trapezoidal rule
		 }
		 printf("%lf",result);
	}
	else if(choice == 10){
		printf("Enter starting x value (x0): ");
		scanf("%lf",&a); // a == x0
		printf("Enter gap (h): ");
		scanf("%lf",&h);
		printf("Enter the count of datas you have (n): ");
		scanf("%d",&n);
		printf("Enter values of these datas.\n");
		
		temp = a;
		for(i = 0; i < n; i++){
			printf("f(%lf): ",temp);
			scanf("%lf",&dataValues[i]);
			temp += h;
		}
		printf("Enter the x value that you want to calculate: ");
		scanf("%lf",&x);
		
		for(i = 1; i < n; i++){
			temp = a;
			entResult = 1;
			for(j = 1; j <= i; j++){
				entResult = entResult*(x - temp)/(h*j);
				temp += h;
			}
			entResult = entResult*forwardDiff(dataValues,i,0);
			result += entResult;
		}
		result = result + dataValues[0];
		printf("Result: %lf",result);
	}
    return 0;
    }




