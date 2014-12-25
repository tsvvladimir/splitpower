/**
 * Created by vladimirtsvetkov on 25/12/14.
 */
public class Sorter {
    public double[] mas;
    public int[] perm;
    public Sorter(double[] inp) {
        this.mas = inp;
        for(int i = 0; i < inp.length; i++) {
            perm[i] = i;
        }
    }
    public double[] perfsort() {
        boolean swapped = true;
        int j = 0;
        double tmp;
        int tmpind;
        while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < mas.length - j; i++) {
                if (mas[i] > mas[i + 1]) {
                    tmp = mas[i];
                    mas[i] = mas[i + 1];
                    mas[i + 1] = tmp;
                    swapped = true;
                    tmpind = perm[i];
                    perm[i] = perm[i + 1];
                    perm[i + 1] = tmpind;
                }
            }
        }
    }
}
