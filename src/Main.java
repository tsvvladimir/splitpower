public class Main {
    public static void main(String[] args) {
        //Matrix m = Matrix.random3diag(4, 4);
        //Matrix m = new Matrix(new double[][] {  {6, 3.5, 0},
        //        {3.5, 5, 2},
        //        {0, 2, 7}});
        Matrix m = new Matrix(new double[][] {  {1, 2, 3},
                                                {2, 1, 2},
                                                {3, 2, 1}});
        //Matrix m = Matrix.identity(5);
        m.show();
        Pow pow = new Pow(m);
        PairMatrix res = pow.p;
        System.out.println("first:\n");
        res.a.show2();
        System.out.println("second:\n");
        res.b.show2();
        System.out.println("пока");
    }
}

//vekovoye: 0 = 1 + 2 * (1/(5-x) + 1/(-0.5-x))