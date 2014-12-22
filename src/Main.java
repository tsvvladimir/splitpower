public class Main {
    public static void main(String[] args) {
        //Matrix m = Matrix.random3diag(4, 4);
        Matrix m = new Matrix(new double[][] {  {6, 3.5, 0},
                                                {3.5, 5, 2},
                                                {0, 2, 7}});
        m.show();
        Pow pow = new Pow(m);
        PairMatrix res = pow.p;
        System.out.println("first:\n");
        res.a.show();
        System.out.println("second:\n");
        res.b.show();
        System.out.println("пока");
    }
}
