/*************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones immutable data type for M-by-N matrices.
 *
 *************************************************************************/

final public class Matrix {
    public final int M;             // number of rows
    public final int N;             // number of columns
    public final double[][] data;   // M-by-N array

    // create M-by-N matrix of 0's
    public Matrix(int M, int N) {
        this.M = M;
        this.N = N;
        data = new double[M][N];
        for(int i = 0; i < M; i++)
            for(int j = 0; j < N; j++)
                data[i][j] = 0.0;
    }


    // create matrix based on 2d array
    public Matrix(double[][] data) {
        M = data.length;
        N = data[0].length;
        this.data = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                this.data[i][j] = data[i][j];
    }

    public Matrix(Double[][] data) {
        M = data.length;
        N = data[0].length;
        this.data = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                this.data[i][j] = data[i][j];
    }

    // copy constructor
    public Matrix(Matrix A) { this(A.data); }

    static public Matrix rowToColumn(Double []a) {
        Matrix m = new Matrix(a.length, 1);
        for (int i = 0; i < a.length; i++)
            m.setElement(i, 0, a[i]);
        return m;
    }



    static public Matrix rowToDiag(Double []a) {
        Matrix m = new Matrix(a.length, a.length);
        for (int i = 0; i < a.length; i++)
            m.setElement(i, i, a[i]);
        return m;
    }

    //changed
    public double GetElement(int i, int j) {
        return data[i][j];
    }

    //changed
    public void setElement(int i, int j, double value) {
        data[i][j] = value;
    }

    // create and return a random M-by-N matrix with values between 0 and 1
    public static Matrix random(int M, int N) {
        Matrix A = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[i][j] = Math.random();
        return A;
    }

    //create and return sym random matrix
    public static Matrix randomsym(int M, int N) {
        Matrix A = Matrix.random(M, N);
        System.out.println("rand before sym");
        A.show();
        Matrix B = new Matrix(A);
        for (int i = 0; i < (M - 1); i ++)
            for(int j = i + 1; j < N; j++)
                B.setElement(i, j, A.GetElement(j, i));
        return B;
    }

    public static Matrix random3diag(int M, int N) {
        Matrix A = Matrix.randomsym(M, N);
        for (int i = 2; i < M; i++)
            for (int j = 0; j <= ( - 2 + i); j++) {
                A.setElement(i, j, 0);
                //System.out.println(i + "--" + j);
            }
        for (int i = 2; i < M; i++)
            for (int j = 0; j <= ( - 2 + i); j++) {
                A.setElement(j, i, 0);
                //System.out.println(i + "--" + j);
            }

        System.out.println("3diag Matrix:\n");
        A.show();
        return A;
    }

    public static Matrix Hilbert(int M, int N){
        Matrix A = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for(int j = 0; j < N; j++) {
                A.data[i][j] = 1 / (double)(i + j + 1);
                //System.out.println("1 / (" + i + "+" + j + "+" + "1 ) = " + (1/(i+j+1))  );
            }
        //System.out.println("\nHilbert Matrix: \n");
        //A.show();
        //System.out.println("\n \n");
        return A;
    }

    // create and return the N-by-N identity matrix
    public static Matrix identity(int N) {
        Matrix I = new Matrix(N, N);
        for (int i = 0; i < N; i++)
            I.data[i][i] = 1;
        return I;
    }

    public static Matrix Null(int M, int N) {
        double [][]a = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                a[i][j] = 0;
        return new Matrix(a);
    }

    public Matrix getBlock(int fromRaw, int toRaw, int fromColumn, int toColumn) {
        if(data == null || fromRaw > toRaw || fromColumn > toColumn || toRaw-1 > data.length || toColumn-1 > data[0].length) {
            System.out.println("error in Matrix.getBlock - can't get such block");
            return null;
        }
        int n = toRaw-fromRaw, m = toColumn-fromColumn;
        double [][]a = new double[n][m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                a[i][j] = data[fromRaw+i][fromColumn+j];
        return new Matrix(a);
    }

    // swap rows i and j
    public void swap(int i, int j) {
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }

    //changed
    public void swapRaws(int i, int j) {
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }
    public void swapColumns(int i, int j) {
        for(int z = 0; z < N; z++) {
            double temp = data[z][i];
            data[z][i] = data[z][j];
            data[z][j] = temp;
        }
    }


    // create and return the transpose of the invoking matrix
    public Matrix transpose() {
        Matrix A = new Matrix(N, M);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }

    //return C = A * digit
    public Matrix muldig(double a) {
        Matrix C = new Matrix(this);
        for(int i = 0; i < C.M; i++)
            for(int j = 0; j < C.N; j++)
                C.setElement(i, j, C.GetElement(i, j) * a);
        return C;
    }

    //return C = A * digit
    public Matrix muldiag(double a) {
        Matrix C = new Matrix(this);
        for(int i = 0; i < C.M; i++)
            C.setElement(i, i, C.GetElement(i, i) * a);
        return C;
    }

    // return C = A + B
    public Matrix plus(Matrix B) {
        Matrix A = new Matrix(this);
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }


    // return C = A - B
    public Matrix minus(Matrix B) {
        Matrix A = new Matrix(this);
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        return C;
    }

    // does A = B exactly?
    public boolean eq(Matrix B) {
        Matrix A = new Matrix(this);
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                if (A.data[i][j] != B.data[i][j]) return false;
        return true;
    }

    // return C = A * B
    public Matrix times(Matrix B) {
        Matrix A = new Matrix(this);

        if (A.N != B.M) throw new RuntimeException("Illegal matrix dimensions." + A.M + " " + A.N + " " + B.M + " " + B.N);
        Matrix C = new Matrix(A.M, B.N);
        for (int i = 0; i < C.M; i++)
            for (int j = 0; j < C.N; j++)
                for (int k = 0; k < A.N; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }

    public Matrix mulOnConstant(double constant) {
        if (data == null) { return null; }
        double [][]a = new double[data.length][data[0].length];
        for (int i = 0; i < data.length; i++)
            for (int j = 0; j < data[i].length; j++)
                a[i][j] = constant * data[i][j];
        return new Matrix(a);
    }


    // return x = A^-1 b, assuming A is square and has full rank
    public Matrix solve(Matrix rhs) {
        if (M != N || rhs.M != N || rhs.N != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this);
        Matrix b = new Matrix(rhs);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < N; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < N; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < N; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < N; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < N; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        // back substitution
        Matrix x = new Matrix(N, 1);
        for (int j = N - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < N; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        return x;

    }

    //changed
    public Matrix unMinus() {
        Matrix a = new Matrix(this);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                a.data[i][j] = -this.data[i][j];
        return a;
    }

    public Matrix getNormed() {
        Matrix a = new Matrix(this);
        Double norm = 0.0;
        for (int i = 0; i < M; i++) {
            norm += a.GetElement(i, 0) * a.GetElement(i, 0);
        }
        norm = Math.sqrt(norm);
        for (int i = 0; i < M; i++) {
            a.setElement(i, 0, a.GetElement(i, 0) / norm);
        }
        if (a.GetElement(0, 0) < 0) {
            //a.muldig(-1.0);
        }
        return a;
    }

    //changed
    // return x = A^-1 , assuming A is square and has full rank
    public Matrix degMin1() {
        double[][] A = new Matrix(this.data).data;
        int n = N;
        int row[] = new int[n];
        int col[] = new int[n];
        double temp[] = new double[n];
        int hold, I_pivot, J_pivot;
        double pivot, abs_pivot;

        if (A[0].length != n) {
            System.out.println("Error in Matrix.invert, inconsistent array sizes.");
        }
        // установиим row и column как вектор изменений.
        for (int k = 0; k < n; k++) {
            row[k] = k;
            col[k] = k;
        }
        // начало главного цикла
        for (int k = 0; k < n; k++) {
            // найдем наибольший элемент для основы
            pivot = A[row[k]][col[k]];
            I_pivot = k;
            J_pivot = k;
            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    abs_pivot = Math.abs(pivot);
                    if (Math.abs(A[row[i]][col[j]]) > abs_pivot) {
                        I_pivot = i;
                        J_pivot = j;
                        pivot = A[row[i]][col[j]];
                    }
                }
            }
            if (Math.abs(pivot) < 1.0E-20) {
                System.out.println("Matrix is singular !");
                return new Matrix(this);
            }
            //перестановка к-ой строки и к-ого столбца с стобцом и строкой, содержащий основной элемент(pivot основу)
            hold = row[k];
            row[k] = row[I_pivot];
            row[I_pivot] = hold;
            hold = col[k];
            col[k] = col[J_pivot];
            col[J_pivot] = hold;
            // k-ую строку с учетом перестановок делим на основной элемент
            A[row[k]][col[k]] = 1.0 / pivot;
            for (int j = 0; j < n; j++) {
                if (j != k) {
                    A[row[k]][col[j]] = A[row[k]][col[j]] * A[row[k]][col[k]];
                }
            }
            // внутренний цикл
            for (int i = 0; i < n; i++) {
                if (k != i) {
                    for (int j = 0; j < n; j++) {
                        if (k != j) {
                            A[row[i]][col[j]] = A[row[i]][col[j]] - A[row[i]][col[k]] *
                                    A[row[k]][col[j]];
                        }
                    }
                    A[row[i]][col[k]] = -A[row[i]][col[k]] * A[row[k]][col[k]];
                }
            }
        }
        // конец главного цикла

        // переставляем назад rows
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                temp[col[i]] = A[row[i]][j];
            }
            for (int i = 0; i < n; i++) {
                A[i][j] = temp[i];
            }
        }
        // переставляем назад columns
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp[row[j]] = A[i][col[j]];
            }
            for (int j = 0; j < n; j++) {
                A[i][j] = temp[j];
            }
        }

        return new Matrix(A);

    }

    //count norm of matrix as sum of its elemnts divided on number of elements
    public double norm1() {
        double result = 0;
        for (int i = 0; i < data.length; i++)
            for (int j = 0; j < data[0].length; j++)
                result += Math.abs(data[i][j]);
        return result/(data.length*data[0].length);
    }
    //counts norm of matrixs as sqrt of sum of elements^2
    public double norm2() {
        double result = 0;
        for (int i = 0; i < data.length; i++)
            for (int j = 0; j < data[0].length; j++)
                result += data[i][j]*data[i][j];
        return Math.sqrt(result);
    }

    // print matrix to standard output
    public void show() {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++)
                System.out.print("(" + data[i][j] + ") ");//f("%f ", data[i][j]);//("%9.10f ", data[i][j]);
            //System.out.printf("%9.10f ", data[i][j]);
            System.out.println();
        }
    }

    // print matrix to standard output
    public void show2() {
        System.out.print("{");
        for (int i = 0; i < M; i++) {
            System.out.print("{");
            for (int j = 0; j < N; j++)
                System.out.printf("%4.4f;", data[i][j]);//f("%f ", data[i][j]);//("%9.10f ", data[i][j]);
            //System.out.printf("%9.10f ", data[i][j]);
            System.out.print("},");
        }
        System.out.println("}");
    }


    // test client
    public static void main(String[] args) {
        double[][] d = { { 1, 2, 3 }, { 4, 5, 6 }, { 9, 1, 3} };
        Matrix D = new Matrix(d);
        D.show();
        System.out.println();

        Matrix A = Matrix.random(5, 5);
        A.show();
        System.out.println();

        A.swap(1, 2);
        A.show();
        System.out.println();

        Matrix B = A.transpose();
        B.show();
        System.out.println();

        Matrix C = Matrix.identity(5);
        C.show();
        System.out.println();

        A.plus(B).show();
        System.out.println();

        B.times(A).show();
        System.out.println();

        // shouldn't be equal since AB != BA in general
        System.out.println(A.times(B).eq(B.times(A)));
        System.out.println();

        Matrix b = Matrix.random(5, 1);
        b.show();
        System.out.println();

        Matrix x = A.solve(b);
        x.show();
        System.out.println();

        A.times(x).show();

    }
}