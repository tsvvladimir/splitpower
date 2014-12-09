public class Main {
    public static void main(String[] args) {
        Matrix m = Matrix.random3diag(4, 4);
        Pow pow = new Pow(m);
        PairMatrix res = pow.p;
        System.out.println("first:\n");
        res.a.show();
        System.out.println("second:\n");
        res.b.show();
    }
}
