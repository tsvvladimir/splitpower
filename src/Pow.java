/**
 * Created by vladimirtsvetkov on 09/12/14.
 */
public class Pow {
    public PairMatrix p;
    public Pow(Matrix T){
        p = dc_eig(T);
    }
    public PairMatrix dc_eig(Matrix T) {
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

            Pow pow1 = new Pow(T1);
            PairMatrix res1 = pow1.p;

            Pow pow2 = new Pow(T2);
            PairMatrix res2 = pow2.p;

            //buildD
            Matrix D = new Matrix(res1.b.M + res2.b.M, res1.b.M + res2.b.M);

            for (int i = 0; i < res1.b.M; i++)
                for(int j = 0; j < res1.b.M; j++)
                    D.setElement(i, j, res1.b.GetElement(i, j));
            for (int i = 0; i < res2.b.M; i++)
                for(int j = 0; j < res2.b.M; j++)
                    D.setElement(i + res1.b.M, j + res1.b.M, res2.b.GetElement(i, j));

            //build u
            Matrix Q1t = res1.a.transpose();
            Matrix Q2t = res2.a.transpose();

            Matrix u = new Matrix(res1.b.M + res2.b.M, 1);
            for (int i = 0; i < res1.b.M; i++)
                u.setElement(i, 0, Q1t.GetElement(i, 0));

            for (int i = 0; i < res2.b.M; i++)
                u.setElement(i + res1.b.M, 0, Q2t.GetElement(i, 0));




            QL.a = T1;
            QL.b = T2;
        }
        return QL;
    }
}
