import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

/*
 * https://cstheory.stackexchange.com/questions/33093/min-hamming-distance-of-a-given-string-from-substrings-of-another-string
 */
class AdaAndNucleobase {

	private static String BASE = "ACTG";
	private static double eps = 1e-3;
	private final static int MAX = 500001;
	private static int[] possibleCandidates = new int[2 * MAX];

	static class Complex {

		public double real;
		public double imaginary;

		public Complex() {
			this.real = 0.0;
			this.imaginary = 0.0;
		}

		public Complex(double r, double i) {
			this.real = r;
			this.imaginary = i;
		}

		public Complex multiply(final Complex x) {
			final Complex copy = new Complex(this.real, this.imaginary);
			copy.real = this.real * x.real - this.imaginary * x.imaginary;
			copy.imaginary = this.imaginary * x.real + this.real * x.imaginary;
			return copy;
		}

		public Complex add(final Complex x) {
			final Complex copy = new Complex(this.real, this.imaginary);
			copy.real += x.real;
			copy.imaginary += x.imaginary;
			return copy;
		}

		public Complex sub(final Complex x) {
			final Complex copy = new Complex(this.real, this.imaginary);
			copy.real -= x.real;
			copy.imaginary -= x.imaginary;
			return copy;
		}

		public double abs() {
			return Math.sqrt(this.real * this.real + this.imaginary * this.imaginary);
		}

		public String toString() {
			return "(" + this.real + "," + this.imaginary + ")";
		}

		public static Complex polar(final double rho, final double theta) {
			return (new Complex(rho * Math.cos(theta), rho * Math.sin(theta)));
		}

		public Complex conjugate() {
			final Complex copy = new Complex(this.real, this.imaginary);
			copy.imaginary = copy.imaginary * -1;
			return copy;

		}

		public void print() {
			final Complex copy = new Complex(this.real, this.imaginary);
			System.out.print("( " + copy.real + "," + copy.imaginary + ") ,");

		}

		public Complex scal(Double dFactor) {
			final Complex copy = new Complex(this.real, this.imaginary);

			copy.real = dFactor * copy.real;
			copy.imaginary = dFactor * copy.imaginary;

			return (copy);
		}

		public Complex Div(Complex cB) {
			Complex div = new Complex();
			double dR, dDen;

			if (Math.abs(cB.real) >= Math.abs(cB.imaginary)) {
				dR = cB.imaginary / cB.real;
				dDen = cB.real + dR * cB.imaginary;
				div.real = (this.real + dR * this.imaginary) / dDen;
				div.imaginary = (this.imaginary - dR * this.real) / dDen;
			} else {
				dR = cB.real / cB.imaginary;
				dDen = cB.imaginary + dR * cB.real;
				div.real = (dR * this.real + this.imaginary) / dDen;
				div.imaginary = (dR * this.imaginary - this.real) / dDen;
			}

			return (div);
		}
	}

	public static Complex[] fft(Complex[] poly, boolean invert) {

		int n = poly.length;

		for (int i = 1, j = 0; i < n; i++) {

			int bit = n >> 1;
			for (; (j & bit) != 0; bit >>= 1) {
				j ^= bit;
			}
			j ^= bit;

			if (i < j) {
				Complex tmp = new Complex();
				tmp = poly[i];
				poly[i] = poly[j];
				poly[j] = tmp;

			}
		}

		for (int len = 2; len <= n; len <<= 1) {

			double ang = 2 * Math.PI / len * (invert ? -1 : 1);
			Complex wlen = new Complex(Math.cos(ang), Math.sin(ang));

			for (int i = 0; i < n; i += len) {
				Complex w = new Complex(1, 0.0);

				for (int j = 0; j < len / 2; j++) {

					Complex u = poly[i + j], v = poly[i + j + len / 2].multiply(w);

					poly[i + j] = u.add(v);
					poly[i + j + len / 2] = u.sub(v);
					w = w.multiply(wlen);
				}

			}

		}

		if (invert) {
			for (int i = 0; i < poly.length; i++)
				poly[i] = poly[i].Div(new Complex(n, 0));
		}

		return poly;

	}

	public static int[] multiply(int[] a, int[] b) {

		int n = 1;
		while (n < a.length + b.length)
			n <<= 1;
		Complex[] fa = new Complex[n];
		Complex[] fb = new Complex[n];

		Arrays.fill(fa, new Complex());
		Arrays.fill(fb, new Complex());

		for (int i = 0; i < a.length; i++) {
			fa[i] = new Complex(a[i], 0.0);
		}

		for (int i = 0; i < b.length; i++) {
			fb[i] = new Complex(b[i], 0.0);
		}

		fa = fft(fa, false);
		fb = fft(fb, false);

		for (int i = 0; i < n; i++) {
			fa[i] = fa[i].multiply(fb[i]);
		}

		fa = fft(fa, true);

		int[] res = new int[n];
		for (int i = 0; i < n; i++) {
			res[i] = (int) Math.round(fa[i].real + eps);
		}

		return res;
	}

	public static void main(String[] args) throws IOException {

		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		int[] ams;
		int[] bss;
		String MS = br.readLine().strip();
		String SS = br.readLine().strip();
		int msLen = 0;
		int ssLen = 0;

		for (int i = 0; i < 4; i++) {

			msLen = MS.length();
			ssLen = SS.length();
			ams = new int[msLen];
			bss = new int[ssLen];
			for (int j = 0; j < msLen; j++) {
				ams[j] = (MS.charAt(j) == BASE.charAt(i)) ? 1 : 0;
			}
			for (int j = ssLen - 1; j >= 0; j--) {
				bss[ssLen - 1 - j] = (SS.charAt(j) == BASE.charAt(i)) ? 1 : 0;
			}
			int[] response = multiply(ams, bss);
			for (int j = 0; j < response.length; j++) {
			}

			for (int j = 0; j < msLen; j++) {
				possibleCandidates[j + ssLen - 1] += response[j + ssLen - 1];
			}
		}
		int maxAns = Integer.MIN_VALUE;
		for (int i = ssLen - 1; i < msLen; i++) {
			maxAns = Math.max(maxAns, possibleCandidates[i]);
		}
		System.out.println(ssLen - maxAns);
	}

}
