public class Gauss {



	public static void main(String[] args) {
		double alpha1, beta1, lambda1 = -1.0, lambda2 = -1.0,rhs1 = 1.3243, rhs2 = 2.4756;
		alpha1 = lambda1; 
		beta1 = rhs1;

		double aks[] = new double [17]; 	// we are creating array as much as equation that is 18 this is for U_k-1
		double bks[] = new double [17];		// we are creating array as much as equation that is 18 this is for U_k
		double cks[] = new double [17];		// we are creating array as much as equation that is 18	this is for U_k+1
		double rhs[] = new double [17];	// we are creating array as much as equation that is 18	this is for right side of equation

		double alphas[] = new double [17];	// array for alphas
		alphas[1] = alpha1;					// first alpha value is equal to lambda
		double betas[] = new double [17];	// array for betas
		betas[1] = beta1;					// first beta value equal to rhs of first equation
		double roots [] = new double [17]; 	// this array roots of the equation 

		aks = calculateAks(aks);			// this calling the method that method identify of coefficient of u_k-1
		bks = calculateBks(bks);			// this calling the method that method identify of coefficient of u_k
		cks = calculateCks(cks);			// this calling the method that method identify of coefficient of u_k+1
		rhs = calculateRhs(rhs);			// this calling the method that identify of rhs of every equation
		
		alphas = calculateAlphas(alphas, calculateBks(bks), calculateCks(cks), calculateAks(aks)); 			//with variables(b, c and a), its find alphas
		betas = calctulateBetas(betas, calculateRhs(rhs), calculateAks(aks), calculateCks(cks), calculateAlphas(alphas, bks, cks, aks)); 		//with variables(rhs, aks, alpha and cks), its find betas
		
		roots [16] = calculateLastRoot(rhs2, lambda2, betas[16], alphas[16]);	//this calling the method that find last root it is u 15
		roots = calculateRoots(roots, alphas, betas);							//this calling the method that find rest of all undefined u
		roots[1] = calculateFirstRoot(lambda1, roots[2], rhs1);					// this calling the method that find first root, it is u 0

		//All roots are found !
		for (int i = 1; i < roots.length; i++) {								//print all roots with values
			System.out.println("U " + (i - 1) + " : " + roots[i]);
		}

	}
	public static double calculateFirstRoot (double lambda1, double root1, double myu1){		//its find first root
		double firstRoot = (lambda1 * root1) + myu1;
		return firstRoot;
	}
	public static double calculateLastRoot (double myu2, double lambda2, double betasLast, double alphaLast){ //its find last root
		double lastRoot = ((myu2) + (lambda2 * betasLast)) / ((1) - (alphaLast * lambda2));
		return lastRoot;
	}
	public static double [] calculateRoots (double arr[], double alphas[], double betas[]){

		int i = 15;
		while (i > 1)
		{
			arr[i] = (alphas[i + 1] * arr[i + 1]) + betas[i + 1];
			i--;
		}
		return arr;
	}
	public static double [] calculateAks (double arr[]){		//this method find coefficient of u_k-1
		for (int i = 1; i < arr.length; i++) {
			arr[i] = Math.exp(-1*(0.3)*i) ;
		}
		return arr;
	}
	public static double [] calculateBks(double arr[]){			//this method find coefficient of u_k
		for (int i = 1; i < arr.length; i++) {
			arr[i] = Math.exp(-1*(0.6)*i);
		}
		return arr;
	}
	public static double [] calculateCks (double arr[]){		//this method find coefficient of u_k+1
		for (int i = 1; i < arr.length; i++) {
			arr[i] = 1;
		}
		return arr;
	}
	public static double [] calculateRhs (double arr[]){		// this method identify of rhs of every equation
		for (int i = 1; i < arr.length; i++) {
			arr[i] = -1*Math.exp(-1*(0.36));
		}
		return arr;
	}
	public static double [] calculateAlphas (double arr[], double bs[], double cs[], double as[]){ 		//this method give alphas value from
		for (int i = 2; i < arr.length; i++) {															// alpha(k+1)= b(k)/(c(k)-a(k)*alpha(k))	
			arr[i] = (bs[i -1]) / ((cs[i - 1]) - ((as[i - 1]) * (arr [i - 1])));
		}
		return arr;
	}
	public static double [] calctulateBetas (double arr[], double rhs[], double aks[], double cks[], double alphas[]){		//this method give alphas value from
		for (int i = 2; i < arr.length; i++) {																					//beta(k+1)= (mu(k)+a(k)*beta(k))/(c(k)-a(k)*alpha(k))
			arr[i] = ((rhs[i - 1]) + (aks[i - 1] * arr[i - 1])) / ((cks[i - 1]) - (aks[i - 1] * alphas[i - 1]));
		}
		return arr;
	}
	
	
}
