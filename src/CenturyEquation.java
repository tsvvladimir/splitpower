import java.util.Arrays;
import java.util.LinkedList;

/**
 * Created by vladimirtsvetkov on 22/12/14.
 */
public class CenturyEquation {
    int k;
    double ro;//, lambda;
    double []u, d;

    int []cases;

    public double [] reverse(double []a) {
        for (int k = 0; k < a.length/2; k++) {
            double temp = a[k];
            a[k] = a[a.length-(1+k)];
            a[a.length-(1+k)] = temp;
        }
        return a;
    }

    public CenturyEquation(int k, double ro, double []u, double []d) { //, double lambda){
        this.k = k;
        this.ro = ro;
        //this.lambda = lambda;
        this.u = u;
        this.d = d;
        cases = new int[this.d.length];
        //Arrays.sort(this.d);
        //this.d = reverse(this.d);

    }

    public CenturyEquation(int k, double ro, Matrix u, Matrix d) {
        this.k = k;
        this.ro = ro;

        this.d = new double[d.N];
        for (int i = 0; i < d.N; i++)
            this.d[i] = d.GetElement(i, i);
        //Arrays.sort(this.d);
        //this.d = reverse(this.d);
        cases = new int[this.d.length];

        this.u = new double[u.M];
        for (int i = 0; i < u.M; i++)
            this.u[i] = u.GetElement(i, 0);

        if (false) {
            System.out.println("d - must be diagonal: ");
            d.show();
            System.out.println("u - must be a column");
            u.show();

            for (int i = 0; i < this.d.length; i++)
                System.out.print(" " + this.d[i]);
            System.out.println();
        }
    }

    //------------for---interval---(d[i+1],d[i])-----------------------
    private double psi_1_sh(double lambda_j, int i) {
        double psi1sh = 0;
        for (int j = 0; j <= i; j++)
            //psi1sh += ro*u[j]/((d[j] - lambda_j)*(d[j] - lambda_j));//
            psi1sh += ro*u[j]*u[j]/((d[j] - lambda_j)*(d[j] - lambda_j));
        //System.out.println("psi1sh: " + psi1sh);
        return psi1sh;
    }

    private double psi_1(double lambda_j, int i) {
        double psi1 = 0;
        for (int j = 0; j <= i; j++)
            //psi1 += ro*u[j]/(d[j] - lambda_j);
            psi1 += ro*u[j]*u[j]/(d[j] - lambda_j);
        //System.out.println("psi1: " + psi1);
        return psi1;
    }


    private double psi_2_sh(double lambda_j, int i) {
        double psi1sh = 0;
        for (int j = i + 1; j < u.length; j++)
        {
            //psi1sh += ro*u[j]/((d[j] - lambda_j)*(d[j] - lambda_j));//
            psi1sh += ro*u[j]*u[j]/((d[j] - lambda_j)*(d[j] - lambda_j));
            //System.out.println("in psi2sh dj: " + d[j] + " lambdaj: " + lambda_j);
        }

       // System.out.println("psi2sh: " + psi1sh);
        return psi1sh;
    }

    private double psi_2(double lambda_j, int i) {
        double psi1 = 0;
        for (int j = i + 1; j < u.length; j++)
            //psi1 += ro*u[j]/(d[j] - lambda_j);//
            psi1 += ro*u[j]*u[j]/(d[j] - lambda_j);
        //System.out.println()
        //System.out.println("psi2: " + psi1);
        return psi1;

    }
    private double interval(int i) {
        double di = d[i], diplus1 = d[i+1];
        double lambda0 = (di + diplus1)/2;

        for (int iter = 0; iter < k; iter++) {
            double psi1sh = psi_1_sh(lambda0, i);
            double c1 = psi1sh * (d[i] - lambda0)*(d[i] - lambda0);
            double c1_kr = psi_1(lambda0, i) - psi1sh * (d[i] - lambda0);


            double psi2sh = psi_2_sh(lambda0, i);
            double c2 = psi2sh * (d[i + 1] - lambda0)* (d[i + 1] - lambda0);
            //System.out.println(psi_2(lambda0, i) + " and " + psi1sh * (d[i + 1] - lambda0));
            double c2_kr = psi_2(lambda0, i) - psi2sh * (d[i + 1] - lambda0);

            //System.out.println("ckr1: " + c1_kr + " ckr2: " + c2_kr);
            double c3 = c1_kr + c2_kr + 1;

            double a = c3;
            //System.out.println("try get " + i + " and length d is " + d.length);
            double b = (-c3) * d[i] - c3 * d[i + 1] - c1 - c2;
            double c = c3 * d[i] * d[i + 1] + c1 * d[i + 1] + c2 * d[i];
            //System.out.println(" a: " + a + " b: " + b + " c: " + c);

            double D = b * b - 4 * a * c;
            if (D < 0) {
                System.out.println("!!! Error: D < 0");
                return 0;
            }
            double x1 = (-b + Math.sqrt(D)) / (2 * a);
            double x2 = (-b - Math.sqrt(D)) / (2 * a);
            double result = 0;
            if (x1 >= d[i + 1] && x1 <= d[i]) {
                result = x1;
            } else if (x2 >= d[i + 1] && x2 <= d[i]) {
                result = x2;
            } else {
                System.out.println("unexpected solution of h(x) = 0; with x1:" + x1 + " and x2:" + x2 + " with di+1:" + d[i + 1] + " and di:" + d[i] );
            }
            lambda0 = result;
        }
        return lambda0;
    }

    //-----------------for---tail-------------------------
    private double f_sh(double lambda) {
        double f = 0;

        double sum = 0;
        for (int i = 0; i < u.length; i++) {
            sum += ((u[i]*u[i])/((d[i]-lambda)*(d[i]-lambda)));
        }

        f += ro*sum;
        return f;
    }

    private double f(double lambda) {
        double f = 1;

        double sum = 0;
        for (int i = 0; i < u.length; i++) {
            sum += (u[i]*u[i]/(d[i]-lambda));
        }

        f += ro*sum;
        return f;
    }

    private double tail(int i, double lambda0) {
        double lambda = lambda0;

        for (int iter = 0; iter < k; iter++) {
            double fsh = f_sh(lambda);

            double c1 =(d[i] - lambda)* (d[i] - lambda) * fsh;
            double c2 = f(lambda) - fsh*(d[i] - lambda);

            lambda = d[i] + c1/c2;//c2/c1;//
        }

        return lambda;
    }

    //finds all eigenvalues for century_equation initialized
    public Double [] count() {
        LinkedList<Double> solutions = new LinkedList<Double>();
        //this.k = 10;//1000;
        boolean addsolflag = true;
        for (int i = 0; i < u.length-1; i++) {
            //System.out.println("try get from u index " + i + " and length u is " + u.length);
            if (u[i] == 0) {
                Double root = d[i];
                if (addsolflag) {
                    System.out.println("add solution with case 1 and root:" + root);
                }
                solutions.add(root);
                cases[i] = 1;
            } else if (d[i] == d[i+1]) {
                Double root = d[i];
                if (addsolflag) {
                    System.out.println("add solution with case 2 and root:" + root);
                }
                solutions.add(root);
                cases[i] = 2;
            } else {
                Double root = interval(i);
                if (addsolflag) {
                    System.out.println("add solution with case 0 and root:" + root);
                }
                solutions.add(root);
                cases[i] = 0;
            }
        }

        if (ro == 0) System.out.println("!!! ERROR - smth goes wrong: ro = 0;");
        if (ro < 0) {
            if (true) {
                System.out.println("start count tail for case ro < 0 for i = " + (u.length-1) + " and lambda0 = " + (d[u.length-1] - 1));
            }
            Double root = tail(u.length-1, d[u.length-1] - 1);
            if (addsolflag) {
                System.out.println("add solution with case ro < 0 and root:" + root);
            }
            solutions.addLast(root);
        } else {
            if (true) {
                System.out.println("start count tail for case ro > 0 for i = " + 0 + " and lambda0 = " + (d[0] + 1));
            }
            Double root = tail(0, d[0] + 1);
            if (addsolflag) {
                System.out.println("add solution with case ro > 0 and root:" + root);
            }
            solutions.addFirst(root);
        }

        Double []result = new Double[solutions.size()];
        for (int i = 0; i < solutions.size(); i++)
            result[i] = solutions.get(i);

        //---------------check solution----------------------
        boolean printFlag = true;//false;//
        if (printFlag) {
            System.out.println("------\ncheck solution of equation: ");
            System.out.println("\tbm: " + ro + ";\n u: ");
            for (int i = 0; i < u.length; i++)
                System.out.print(" " + u[i]);
            System.out.println("\n\tD:");
            for (int i = 0; i < d.length; i++)
                System.out.print(" " + d[i]);
            System.out.println("\nROOTS:");
            for (int i = 0; i < result.length; i++)
                System.out.print(", " + result[i]);
            System.out.println();

            for (int j = 0; j < result.length; j++) {
                double sum = 0;
                for (int i = 0; i < d.length; i++) {
                    sum += (u[i]*u[i] / (d[i] - result[j]));//result[result.length-i-1]);
                }
                double mustBeNol = 1 + ro * sum;
                System.out.println("must be 0: " + mustBeNol);
            }
        }
        //-----------------------------------------------------

        return result;
    }
}