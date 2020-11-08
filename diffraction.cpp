#include <thread>
#include <Windows.h>
#include "RectDiffraction.h"
#include "CircMultiRectDiffraction.h"
#include "LineMultiRectDiffraction.h"
#include "CircDiffraction.h"

LinearAlgebra la;
int threads = 32;
namespace Program {
	int mainMenu();
	void start(matrix& intensity, double a0, double b0, int dima0, int dimb0, void (*initialize)(), complex<double>(*integral)(vector& x2), void (*dispose)());
	void assign(int* task, int** pstn, vector* x2, double a0, double b0, int dima0, int dimb0);
	void compute(matrix* intensity, vector* x2, int* pstn, int length, int* pcount, double a0, double b0, int dima0, int dimb0, complex<double>(*integral)(vector& x2));
	void disposeAllResults();
}
namespace Rect {
	RectDiffraction rd;
	int set();
	void start();
	void initialize();
	complex<double> integral(vector& x2);
	void dispose();
}
namespace CircMultiRect {
	CircMultiRectDiffraction cmrd;
	int set();
	void start();
	void initialize();
	complex<double> integral(vector& x2);
	void dispose();
}
namespace LineMultiRect {
	LineMultiRectDiffraction lmrd;
	int set();
	void start();
	void initialize();
	complex<double> integral(vector& x2);
	void dispose();
}
namespace Circ {
	CircDiffraction cd;
	int set();
	void start();
	void initialize();
	complex<double> integral(vector& x2);
	void dispose();
}

int main() {
	system("title Diffraction");
	Program::mainMenu();
	Program::disposeAllResults();
	return 0;
}

namespace Program {
	int mainMenu() {
		while (true) {
			system("cls");
			cout << "Please select a diffraction type:" << endl;
			cout << "0.Exit the program" << endl;
			cout << "1.Single rectangular aperture diffraction" << endl;
			cout << "2.Circular multi-rectangular aperture diffraction" << endl;
			cout << "3.Linear multi-rectangular aperture diffraction" << endl;
			cout << "4.Single elliptical aperture diffraction" << endl;
			int type = -1;
			cin >> type;
			switch (type) {
			case 1:
				if (Rect::set() == 0) {
					Rect::start();
				}
				break;
			case 2:
				if (CircMultiRect::set() == 0) {
					CircMultiRect::start();
				}
				break;
			case 3:
				if (LineMultiRect::set() == 0) {
					LineMultiRect::start();
				}
				break;
			case 4:
				if (Circ::set() == 0) {
					Circ::start();
				}
				break;
			default:
				return 0;
			}
		}
	}

	void start(matrix& intensity, double a0, double b0, int dima0, int dimb0, void (*initialize)(), complex<double>(*integral)(vector& x2), void (*dispose)()) {
		try {
			cout << endl << "Initializing..." << endl;
			la.mat.initialize(intensity, dimb0, dima0);
			int* task = new int[threads];
			int** pstn = new int* [threads];
			for (int i = 0; i < threads; i++) {
				pstn[i] = new int[2];
			}
			vector* x2 = new vector[threads];
			assign(task, pstn, x2, a0, b0, dima0, dimb0);
			initialize();

			int* count = new int[threads];
			for (int i = 0; i < threads; i++) {
				count[i] = 0;
				thread t(compute, &intensity, x2 + i, pstn[i], task[i], count + i, a0, b0, dima0, dimb0, integral);
				t.detach();
			}
			int N = dima0 * dimb0;
			int total = 0;
			do {
				total = 0;
				for (int i = 0; i < threads; i++) {
					total += count[i];
				}
				float process = (float)total / (float)N * 100;
				system("cls");
				cout << "Computing..." << endl;
				cout << process << '%';
				Sleep(1000);
			} while (total < N);
			la.mat.fprint(intensity);

			delete[]x2;
			for (int i = 0; i < threads; i++) {
				delete[]pstn[i];
			}
			delete[]pstn;
			delete[]task;
			dispose();
			la.mat.dispose(intensity);

			cout << endl << "Completed." << endl;
			system("pause");
		}
		catch (bad_alloc) {
			cout << "Out of memory!" << endl;
			disposeAllResults();
			system("pause");
			exit(-1);
		}
	}

	void assign(int* task, int** pstn, vector* x2, double a0, double b0, int dima0, int dimb0) {
		double dx0 = a0 / (double)dima0;
		double dy0 = b0 / (double)dimb0;
		int N = dima0 * dimb0;
		int length = N / threads;
		for (int i = 0; i < threads - 1; i++) {
			task[i] = length;
		}
		task[threads - 1] = length + N % threads;
		double x0 = -a0 / 2 + dx0 / 2;
		double y0 = b0 / 2 - dy0 / 2;
		int taskCount = 0;
		for (int i = 0; i < threads; i++) {
			x2[i] = NULL_VEC;
			la.vec.initialize(x2[i], 3);
			pstn[i][0] = taskCount / dima0;
			pstn[i][1] = taskCount % dima0;
			x2[i].value[0] = x0 + pstn[i][1] * dx0;
			x2[i].value[1] = y0 - pstn[i][0] * dy0;
			taskCount += task[i];
		}
	}

	void compute(matrix* intensity, vector* x2, int* pstn, int length, int* pcount, double a0, double b0, int dima0, int dimb0, complex<double>(*integral)(vector& x2)) {
		const double dx0 = a0 / (double)dima0;
		const double dy0 = b0 / (double)dimb0;
		const double x0 = -a0 / 2 + dx0 / 2;
		for (int i = pstn[0]; i < dimb0; i++) {
			for (int j = pstn[1]; j < dima0; j++) {
				complex<double>A = integral(*x2);
				double normA = norm(A);
				intensity->value[i][j] = normA;
				(*pcount)++;
				if (*pcount == length) {
					break;
				}
				x2->value[0] += dx0;
			}
			if (*pcount == length) {
				break;
			}
			x2->value[0] = x0;
			x2->value[1] -= dy0;
		}
		la.vec.dispose(*x2);
	}

	void disposeAllResults() {
		la.mat.dispose(Rect::rd.fk.intensity);
		la.mat.dispose(CircMultiRect::cmrd.fk.intensity);
		la.mat.dispose(LineMultiRect::lmrd.fk.intensity);
		la.mat.dispose(Circ::cd.fk.intensity);
	}
}

namespace Rect {
	int set() {
		double A0 = 0, s0 = 0, r0 = 0, a = 0, b = 0, a0 = 0, b0 = 0, lambda = 0;
		int dima = 0, dimb = 0, dima0 = 0, dimb0 = 0;
		while (true) {
			system("cls");
			cout << "Single rectangular aperture diffraction:" << endl;
			cout << "Please follow the instructions to enter each parameter:" << endl;
			cout << "(System of unit: SI)" << endl;
			do {
				cin.clear();
				cin.ignore();
				cout << "Intensity factor:" << endl << "A0 = ";
			} while (!(cin >> A0) || A0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "Wavelength of light:" << endl << "λ = ";
			} while (!(cin >> lambda) || lambda <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between the point light source and the aperture:" << endl << "s0 = ";
			} while (!(cin >> s0) || s0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between the aperture and the screen:" << endl << "r0 = ";
			} while (!(cin >> r0) || r0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the aperture:" << endl << "a = ";
			} while (!(cin >> a) || a <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the aperture:" << endl << "b = ";
			} while (!(cin >> b) || b <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the aperture width:" << endl << "dima = ";
			} while (!(cin >> dima) || dima <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the aperture height:" << endl << "dimb = ";
			} while (!(cin >> dimb) || dimb <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the screen:" << endl << "a0 = ";
			} while (!(cin >> a0) || a0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the screen:" << endl << "b0 = ";
			} while (!(cin >> b0) || b0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen width:" << endl << "dima0 = ";
			} while (!(cin >> dima0) || dima0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen height:" << endl << "dimb0 = ";
			} while (!(cin >> dimb0) || dimb0 <= 0);
			system("cls");
			cout << "Please check the diffraction type and parameters." << endl;
			cout << "Diffraction type: Single rectangular aperture diffraction" << endl;
			cout << "A0    = " << A0 << endl;
			cout << "λ    = " << lambda << endl;
			cout << "s0    = " << s0 << endl;
			cout << "r0    = " << r0 << endl;
			cout << "a     = " << a << endl;
			cout << "b     = " << b << endl;
			cout << "dima  = " << dima << endl;
			cout << "dimb  = " << dimb << endl;
			cout << "a0    = " << a0 << endl;
			cout << "b0    = " << b0 << endl;
			cout << "dima0 = " << dima0 << endl;
			cout << "dimb0 = " << dimb0 << endl;
			cout << "Please enter 0 if everything is OK, and start the computation." << endl;
			cout << "Please enter 1 if you need to modify the parameters." << endl;
			cout << "Please enter 2 if you need to modify the diffraction type." << endl;
			int check = -1;
			do {
				cin.clear();
				cin.ignore();
			} while (!(cin >> check));
			switch (check) {
			case 0:
				rd.fk.A0 = A0;
				rd.fk.lambda = lambda;
				rd.fk.s0 = s0;
				rd.fk.r0 = r0;
				rd.a = a;
				rd.b = b;
				rd.dima = dima;
				rd.dimb = dimb;
				rd.a0 = a0;
				rd.b0 = b0;
				rd.dima0 = dima0;
				rd.dimb0 = dimb0;
				return 0;
			case 1:
				break;
			case 2:
				return 1;
			default:
				return -1;
			}
		}
	}

	void start() {
		Program::start(rd.fk.intensity, rd.a0, rd.b0, rd.dima0, rd.dimb0, initialize, integral, dispose);
	}

	void initialize() {
		rd.initialize();
	}

	complex<double> integral(vector& x2) {
		return rd.integral(x2);
	}

	void dispose() {
		rd.dispose();
	}
}

namespace CircMultiRect {
	int set() {
		double A0 = 0, s0 = 0, r0 = 0, d = 0, a = 0, b = 0, a0 = 0, b0 = 0, lambda = 0;
		int apCount = 0, dima = 0, dimb = 0, dima0 = 0, dimb0 = 0;
		while (true) {
			system("cls");
			cout << "Circular multi-rectangular aperture diffraction:" << endl;
			cout << "Please follow the instructions to enter each parameter:" << endl;
			cout << "(System of unit: SI)" << endl;
			do {
				cin.clear();
				cin.ignore();
				cout << "Intensity factor:" << endl << "A0 = ";
			} while (!(cin >> A0) || A0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "Wavelength of light:" << endl << "λ = ";
			} while (!(cin >> lambda) || lambda <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between the point light source and apertures' center:" << endl << "s0 = ";
			} while (!(cin >> s0) || s0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between apertures' center and the screen:" << endl << "r0 = ";
			} while (!(cin >> r0) || r0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "Total number of apertures:" << endl << "N = ";
			} while (!(cin >> apCount) || apCount <= 1);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between each aperture and apertures' center:" << endl << "d = ";
			} while (!(cin >> d) || d <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the aperture:" << endl << "a = ";
			} while (!(cin >> a) || a <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the aperture:" << endl << "b = ";
			} while (!(cin >> b) || b <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the aperture width:" << endl << "dima = ";
			} while (!(cin >> dima) || dima <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the aperture height:" << endl << "dimb = ";
			} while (!(cin >> dimb) || dimb <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the screen:" << endl << "a0 = ";
			} while (!(cin >> a0) || a0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the screen:" << endl << "b0 = ";
			} while (!(cin >> b0) || b0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen width:" << endl << "dima0 = ";
			} while (!(cin >> dima0) || dima0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen height:" << endl << "dimb0 = ";
			} while (!(cin >> dimb0) || dimb0 <= 0);
			system("cls");
			cout << "Please check the diffraction type and parameters." << endl;
			cout << "Diffraction type: Circular multi-rectangular aperture diffraction" << endl;
			cout << "A0    = " << A0 << endl;
			cout << "λ    = " << lambda << endl;
			cout << "s0    = " << s0 << endl;
			cout << "r0    = " << r0 << endl;
			cout << "N     = " << apCount << endl;
			cout << "d     = " << d << endl;
			cout << "a     = " << a << endl;
			cout << "b     = " << b << endl;
			cout << "dima  = " << dima << endl;
			cout << "dimb  = " << dimb << endl;
			cout << "a0    = " << a0 << endl;
			cout << "b0    = " << b0 << endl;
			cout << "dima0 = " << dima0 << endl;
			cout << "dimb0 = " << dimb0 << endl;
			cout << "Please enter 0 if everything is OK, and start the computation." << endl;
			cout << "Please enter 1 if you need to modify the parameters." << endl;
			cout << "Please enter 2 if you need to modify the diffraction type." << endl;
			int check = -1;
			do {
				cin.clear();
				cin.ignore();
			} while (!(cin >> check));
			switch (check) {
			case 0:
				cmrd.fk.A0 = A0;
				cmrd.fk.lambda = lambda;
				cmrd.fk.s0 = s0;
				cmrd.fk.r0 = r0;
				cmrd.apCount = apCount;
				cmrd.d = d;
				cmrd.a = a;
				cmrd.b = b;
				cmrd.dima = dima;
				cmrd.dimb = dimb;
				cmrd.a0 = a0;
				cmrd.b0 = b0;
				cmrd.dima0 = dima0;
				cmrd.dimb0 = dimb0;
				return 0;
			case 1:
				break;
			case 2:
				return 1;
			default:
				return -1;
			}
		}
	}

	void start() {
		Program::start(cmrd.fk.intensity, cmrd.a0, cmrd.b0, cmrd.dima0, cmrd.dimb0, initialize, integral, dispose);
	}

	void initialize() {
		cmrd.initialize();
	}

	complex<double> integral(vector& x2) {
		return cmrd.integral(x2);
	}

	void dispose() {
		cmrd.dispose();
	}
}

namespace LineMultiRect {
	int set() {
		double A0 = 0, s0 = 0, r0 = 0, d = 0, a = 0, b = 0, a0 = 0, b0 = 0, lambda = 0;
		int apCount = 0, dima = 0, dimb = 0, dima0 = 0, dimb0 = 0;
		while (true) {
			system("cls");
			cout << "Linear multi-rectangular aperture diffraction:" << endl;
			cout << "Please follow the instructions to enter each parameter:" << endl;
			cout << "(System of unit: SI)" << endl;
			do {
				cin.clear();
				cin.ignore();
				cout << "Intensity factor:" << endl << "A0 = ";
			} while (!(cin >> A0) || A0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "Wavelength of light:" << endl << "λ = ";
			} while (!(cin >> lambda) || lambda <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between the point light source and apertures' center:" << endl << "s0 = ";
			} while (!(cin >> s0) || s0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between apertures' center and the screen:" << endl << "r0 = ";
			} while (!(cin >> r0) || r0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "Total number of apertures:" << endl << "N = ";
			} while (!(cin >> apCount) || apCount <= 1);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between adjacent apertures:" << endl << "d = ";
			} while (!(cin >> d) || d <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the aperture:" << endl << "a = ";
			} while (!(cin >> a) || a <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the aperture:" << endl << "b = ";
			} while (!(cin >> b) || b <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the aperture width:" << endl << "dima = ";
			} while (!(cin >> dima) || dima <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the aperture height:" << endl << "dimb = ";
			} while (!(cin >> dimb) || dimb <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the screen:" << endl << "a0 = ";
			} while (!(cin >> a0) || a0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the screen:" << endl << "b0 = ";
			} while (!(cin >> b0) || b0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen width:" << endl << "dima0 = ";
			} while (!(cin >> dima0) || dima0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen height:" << endl << "dimb0 = ";
			} while (!(cin >> dimb0) || dimb0 <= 0);
			system("cls");
			cout << "Please check the diffraction type and parameters." << endl;
			cout << "Diffraction type: Linear multi-rectangular aperture diffraction" << endl;
			cout << "A0    = " << A0 << endl;
			cout << "λ    = " << lambda << endl;
			cout << "s0    = " << s0 << endl;
			cout << "r0    = " << r0 << endl;
			cout << "N     = " << apCount << endl;
			cout << "d     = " << d << endl;
			cout << "a     = " << a << endl;
			cout << "b     = " << b << endl;
			cout << "dima  = " << dima << endl;
			cout << "dimb  = " << dimb << endl;
			cout << "a0    = " << a0 << endl;
			cout << "b0    = " << b0 << endl;
			cout << "dima0 = " << dima0 << endl;
			cout << "dimb0 = " << dimb0 << endl;
			cout << "Please enter 0 if everything is OK, and start the computation." << endl;
			cout << "Please enter 1 if you need to modify the parameters." << endl;
			cout << "Please enter 2 if you need to modify the diffraction type." << endl;
			int check = -1;
			do {
				cin.clear();
				cin.ignore();
			} while (!(cin >> check));
			switch (check) {
			case 0:
				lmrd.fk.A0 = A0;
				lmrd.fk.lambda = lambda;
				lmrd.fk.s0 = s0;
				lmrd.fk.r0 = r0;
				lmrd.apCount = apCount;
				lmrd.d = d;
				lmrd.a = a;
				lmrd.b = b;
				lmrd.dima = dima;
				lmrd.dimb = dimb;
				lmrd.a0 = a0;
				lmrd.b0 = b0;
				lmrd.dima0 = dima0;
				lmrd.dimb0 = dimb0;
				return 0;
			case 1:
				break;
			case 2:
				return 1;
			default:
				return -1;
			}
		}
	}

	void start() {
		Program::start(lmrd.fk.intensity, lmrd.a0, lmrd.b0, lmrd.dima0, lmrd.dimb0, initialize, integral, dispose);
	}

	void initialize() {
		lmrd.initialize();
	}

	complex<double> integral(vector& x2) {
		return lmrd.integral(x2);
	}

	void dispose() {
		lmrd.dispose();
	}
}

namespace Circ {
	int set() {
		double A0 = 0, s0 = 0, r0 = 0, a = 0, b = 0, a0 = 0, b0 = 0, lambda = 0;
		int dimr = 0, dimt = 0, dima0 = 0, dimb0 = 0;
		while (true) {
			system("cls");
			cout << "Single elliptical aperture diffraction:" << endl;
			cout << "Please follow the instructions to enter each parameter:" << endl;
			cout << "(System of unit: SI)" << endl;
			do {
				cin.clear();
				cin.ignore();
				cout << "Intensity factor:" << endl << "A0 = ";
			} while (!(cin >> A0) || A0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "Wavelength of light:" << endl << "λ = ";
			} while (!(cin >> lambda) || lambda <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between the point light source and the aperture:" << endl << "s0 = ";
			} while (!(cin >> s0) || s0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The distance between the aperture and the screen:" << endl << "r0 = ";
			} while (!(cin >> r0) || r0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The major axis of the aperture:" << endl << "a = ";
			} while (!(cin >> a) || a <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The minor axis of the aperture:" << endl << "b = ";
			} while (!(cin >> b) || b <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of radial direction:" << endl << "dimr = ";
			} while (!(cin >> dimr) || dimr <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of angular direction:" << endl << "dimt = ";
			} while (!(cin >> dimt) || dimt <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The width of the screen:" << endl << "a0 = ";
			} while (!(cin >> a0) || a0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The height of the screen:" << endl << "b0 = ";
			} while (!(cin >> b0) || b0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen width:" << endl << "dima0 = ";
			} while (!(cin >> dima0) || dima0 <= 0);
			do {
				cin.clear();
				cin.ignore();
				cout << "The dimension of the screen height:" << endl << "dimb0 = ";
			} while (!(cin >> dimb0) || dimb0 <= 0);
			system("cls");
			cout << "Please check the diffraction type and parameters." << endl;
			cout << "Diffraction type: Single elliptical aperture diffraction" << endl;
			cout << "A0    = " << A0 << endl;
			cout << "λ    = " << lambda << endl;
			cout << "s0    = " << s0 << endl;
			cout << "r0    = " << r0 << endl;
			cout << "a     = " << a << endl;
			cout << "b     = " << b << endl;
			cout << "dimr  = " << dimr << endl;
			cout << "dimt  = " << dimt << endl;
			cout << "a0    = " << a0 << endl;
			cout << "b0    = " << b0 << endl;
			cout << "dima0 = " << dima0 << endl;
			cout << "dimb0 = " << dimb0 << endl;
			cout << "Please enter 0 if everything is OK, and start the computation." << endl;
			cout << "Please enter 1 if you need to modify the parameters." << endl;
			cout << "Please enter 2 if you need to modify the diffraction type." << endl;
			int check = -1;
			do {
				cin.clear();
				cin.ignore();
			} while (!(cin >> check));
			switch (check) {
			case 0:
				cd.fk.A0 = A0;
				cd.fk.lambda = lambda;
				cd.fk.s0 = s0;
				cd.fk.r0 = r0;
				cd.a = a;
				cd.b = b;
				cd.dimr = dimr;
				cd.dimt = dimt;
				cd.a0 = a0;
				cd.b0 = b0;
				cd.dima0 = dima0;
				cd.dimb0 = dimb0;
				return 0;
			case 1:
				break;
			case 2:
				return 1;
			default:
				return -1;
			}
		}
	}

	void start() {
		Program::start(cd.fk.intensity, cd.a0, cd.b0, cd.dima0, cd.dimb0, initialize, integral, dispose);
	}

	void initialize() {
		cd.initialize();
	}

	complex<double> integral(vector& x2) {
		return cd.integral(x2);
	}

	void dispose() {
		cd.dispose();
	}
}