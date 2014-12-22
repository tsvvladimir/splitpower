import java.util.HashMap;

/**
 * Created by vladimirtsvetkov on 09/12/14.
 */
public class Pow {
    public PairMatrix p;


    public Pow(Matrix T){
        p = dc_eig(T, 0.0001, 5);
    }

    public PairMatrix dc_eig(Matrix T, double eps, int k) {
        PairMatrix QL = new PairMatrix(Matrix.identity(1), Matrix.identity(1));
        if (T.M == 1) {
            QL.b = T;
        } else {
            int div = T.M / 2;
            System.out.println("will be div:" + div);
            Matrix T1 = T.getBlock(0, div, 0, div);
            Matrix T2 = T.getBlock(div, T.M, div, T.M);
            double bm = T.GetElement(div - 1, div);
            T1.setElement(T1.M - 1, T1.M - 1, T1.GetElement(T1.M - 1, T1.M - 1) - bm);
            T2.setElement(0, 0, T2.GetElement(0, 0) - bm);
            System.out.println(bm);

            Pow pow1 = new Pow(T1); //first call
            PairMatrix res1 = pow1.p;

            Pow pow2 = new Pow(T2); //second call
            PairMatrix res2 = pow2.p;

            //buildD
            Matrix D = new Matrix(res1.b.M + res2.b.M, res1.b.M + res2.b.M);

            for (int i = 0; i < res1.b.M; i++)
                for(int j = 0; j < res1.b.M; j++)
                    D.setElement(i, j, res1.b.GetElement(i, j));
            for (int i = 0; i < res2.b.M; i++)
                for(int j = 0; j < res2.b.M; j++)
                    D.setElement(i + res1.b.M, j + res1.b.M, res2.b.GetElement(i, j));


            for(int i = 0; i < (res1.b.M + res2.b.M); i++) {
                for (int j = 0; j < (res1.b.M + res2.b.M); j++) {
                    if ((i != j) && ((Math.abs(D.GetElement(i, i) - D.GetElement(j, j))) < eps)) {
                        D.setElement(j, j, D.GetElement(i, i));
                    }
                }
            }

            //build u
            Matrix Q1t = res1.a.transpose();
            Matrix Q2t = res2.a.transpose();

            Matrix u = new Matrix(res1.b.M + res2.b.M, 1);
            for (int i = 0; i < res1.b.M; i++) {
                if (Q1t.GetElement(i, 0) > eps) {
                    u.setElement(i, 0, Q1t.GetElement(i, 0));
                }
            }


            for (int i = 0; i < res2.b.M; i++) {
                if (Q2t.GetElement(i, 0) > eps) {
                    u.setElement(i + res1.b.M, 0, Q2t.GetElement(i, 0));
                }
            }


            //унулить маленькие числа

            //to vekovoe - u matrix and D diagonal

            CenturyEquation cq = new CenturyEquation(k, bm, u, D);
            //have arr []roots
            Double[] roots = cq.count(); //null;
            for(int i = 0 ; i < roots.length; i++)
                System.out.print(", " + roots[i]);
            System.out.println();
            //levner
            //will recieve new w
            HashMap<Double, Integer> getkdij = new HashMap<Double, Integer>();
            for(int i = 0; i < res1.b.M + res2.b.M; i++) {
                if (!getkdij.containsKey(D.GetElement(i, i))) {
                    getkdij.put(D.GetElement(i, i), 1);
                } else {
                    getkdij.put(D.GetElement(i, i), getkdij.get(D.GetElement(i, i)));
                }
            }


            HashMap<Double, Integer> getlj = new HashMap<Double, Integer>();
            for(Double key: getkdij.keySet()) {
                int lol = 0;
                for(double num : roots) {
                    if (key == num){
                        lol++;
                    }
                }
                getlj.put(key, lol);
            }
            Matrix eigenvector = Matrix.Hilbert(1, 1);
            //eigenvector = null;
            double[][] ismQ1 = new double[res1.b.M + res2.b.M][res1.b.M + res2.b.M];
            for (int z = 0; z < res1.b.M + res2.b.M; z++) {
                Double[] w = new Double[res1.b.M + res2.b.M];
                for(int i = 0; i < res1.b.M + res2.b.M; i++) {
                    System.out.println("@@@ " + D.GetElement(i, i) + ", " + getkdij.keySet() + "###");
                    System.out.println("@@" + getkdij.get(D.GetElement(i, i)) + "##");
                    if (getkdij.get(D.GetElement(i, i)) == 1) {
                        System.out.println("ok;");
                        double nom = 1.0;
                        for(int j = 0; j < res1.b.M + res2.b.M; j++) {
                            System.out.println(roots.length);
                            nom *= roots[i] - D.GetElement(i, i);
                        }
                        double denom = 1.0;
                        for(int j = 0; (j < res1.b.M + res2.b.M) ; j++) {
                            if (j != i) {
                                denom *= (D.GetElement(j, j) - D.GetElement(i, i));
                            }
                        }
                        w[i] = Math.sqrt(nom/denom);
                    } else if (getlj.get(roots[i]) > (getkdij.get(roots[i]) - 1)) {
                        w[i] = 0.0;
                    } else if (getlj.get(roots[i]) == (getkdij.get(roots[i]) - 1)) {
                        double nom = 1.0;
                        double denom = 1.0;
                        for(int j = 0; (j < res1.b.M + res2.b.M); j++) {
                            if (roots[j] != D.GetElement(j, j)) {
                                nom *= roots[j] - D.GetElement(j, j);
                            }
                        }
                        for(int j = 0; (j < res1.b.M + res2.b.M); j++) {
                            if (D.GetElement(j, j) != D.GetElement(i, i)) {
                                denom *= getkdij.get(D.GetElement(j, j)) * (D.GetElement(j, j) - D.GetElement(i, i));
                            }
                        }
                    } else {
                        System.out.println("special case");
                    }
                }

                eigenvector = D.minus(Matrix.identity(D.M).muldig(roots[z])).degMin1().times(Matrix.rowToColumn(w));
                for(int y = 0; y < res1.b.M + res2.b.M; y++) {
                    ismQ1[y][z] = eigenvector.GetElement(y, 0);
                }


            }

            //set blocks

            Matrix LLL = Matrix.rowToDiag(roots);

            Matrix Q111 = new Matrix(ismQ1);

            //
            int mm = Q1t.M + Q2t.M, nn = Q1t.N + Q2t.N;
            Matrix Qfinal = new Matrix(Q1t.M + Q2t.M, Q1t.N + Q2t.N);
            for (int z = 0; z < Q1t.M; z++) {
                for (int zz = 0; zz < Q1t.N; zz++)
                    Qfinal.setElement(z, zz, Q1t.GetElement(z, zz));
            }
            for (int z = 0; z < Q2t.M; z++) {
                for (int zz = 0; zz < Q2t.N; zz++)
                    Qfinal.setElement(z + Q1t.M, zz+ Q1t.N, Q2t.GetElement(z, zz));
            }


            QL.a = Qfinal;//Q111;
            QL.b = LLL;
        }
        return QL;
    }
}
