using System;

namespace KURS
{
    class Program
    {
        static double eps = 1e-3;
        static double E = Math.E;

        static void Main()
        {
            double x0 = 0, y0 = 1, x1 = 1, y1 = E, h = 0.2;

            int arr_size = Convert.ToInt32(((x1 - x0) / h) + 1);
            double[,] points = new double[arr_size, 2];

            Console.WriteLine("Исходныеданные:\ny\" = (y' + 2*e^x + y) / 4\ny({0}) = {1:f3}\ny({2}) = {3:f3}\n", x0, y0, x1, y1);
            Console.WriteLine("Находим первую производную в точке x = {0} методом стрельб", x0);
            double firstDerivative = Shooting_Method(x0, y0, x1, y1, h);
            Console.WriteLine("y' = {0:f3}\n", firstDerivative);


            Console.WriteLine("Вычисление значений функий в точках x при помощи метода Рунге-Кутта II порядка с усреднением по времени:");
            for (int i = 0; i < arr_size; i++)
            {
                points[i, 0] = h * i;
                points[i, 1] = Runge_Kutta_Method(x0, x0 + h * i, h, y0, firstDerivative);
                Console.WriteLine("x = {0:f1} y = {1:f3}", points[i, 0], points[i, 1]);
            }

            Console.WriteLine("\nВычисление значений функий в точках x при помощи Кубического сплайна:");
            for (double m = x0; m < x1 + 0.01; m += 0.01)
            {
                Console.WriteLine("x = {0:f2} y = {1:f3}", m, Cubic_Spline(m, arr_size, points));
            }

            Console.WriteLine("\nзначение интеграла:");
            Console.WriteLine(Simpsons_Rule(points, arr_size - 1, firstDerivative, 2));
        }

        // - реализует метод стрельб, который позволяет найти первую производную функции.
        static public double Shooting_Method(double x0, double x1, double y0, double y1, double h)
        {
            double lowerD = 0, upperD = 1;
            double resultA = 1, resultB = 1;

            while (resultA * resultB> 0)
            {
                resultA = Runge_Kutta_Method(x0, x1, h, y0, lowerD) - y1;
                resultB = Runge_Kutta_Method(x0, x1, h, y0, upperD) - y1;

                lowerD -= h;
                upperD += h;
            }

            double g = 0;
            // методом половинного деления ищем первую производную
            while (Math.Abs(upperD - lowerD) > eps)
            {
                g = (lowerD + upperD) / 2;
                if ((Runge_Kutta_Method(x0, x1, h, y0, lowerD) - y1) * (Runge_Kutta_Method(x0, x1, h, y0, g) - y1) < 0)
                {
                    upperD = g;
                }
                else if ((Runge_Kutta_Method(x0, x1, h, y0, g) - y1) * (Runge_Kutta_Method(x0, x1, h, y0, upperD) - y1) < 0)
                {
                    lowerD = g;
                }
            }

            return (lowerD + upperD) / 2;
        }

        // - находит значение функции в точке x1 при помощи метода Рунге-Кутта 2 порядка с усреднением по времени.
        static public double Runge_Kutta_Method(double x0, double x1, double h, double y, double firstDerivative)
        {
            double yDelta, derivativeDelta;
            while (x0 < x1)
            {
                yDelta = y + (h / 2) * firstDerivative; // 1.2
                derivativeDelta = firstDerivative + (h / 2) * Find_Second_Derivative(x0, y, firstDerivative); // 2.2
                y += h * derivativeDelta; // 1.1
                firstDerivative += h * Find_Second_Derivative(x0 + h / 2, yDelta, derivativeDelta); // 2.1
                x0 += h;
            }
            return y;
        }

        // - метод решает уравнение относительно старшей производной и возвращает её значение (второй производной).
        static public double Find_Second_Derivative(double x, double y, double firstDerivative)
        {
            double lowerD = 1, upperD = 0;
            double resultA = 1, resultB = 1;
            while (resultA * resultB> 0)
            {
                resultA = Function(x, y, firstDerivative, lowerD--);
                resultB = Function(x, y, firstDerivative, upperD++);
            }

            doublec = 0;
            // находим вторую производную при помощи МПД
            while (Math.Abs(upperD - lowerD) > eps)
            {
                c = (lowerD + upperD) / 2;
                if (Function(x, y, firstDerivative, lowerD) * Function(x, y, firstDerivative, c) < 0)
                {
                    upperD = c;
                }
                else if (Function(x, y, firstDerivative, c) * Function(x, y, firstDerivative, upperD) < 0)
                {
                    lowerD = c;
                }
                else
                {
                    break;
                }
            }
            return (lowerD + upperD) / 2;
        }

        // - хранит исходное уравнение. Возвращает значение уравнения в точке x и её значениях y, y', y".
        static public double Function(double x, double y, double firstDerivative, double secondDerivative)
        {
            return 4 * secondDerivative - firstDerivative - 2 * Math.Exp(x) - y;
        }

        // - находит длину i-ого интервала.
        static public double Get_Interval(int i, double[,] points)
        {
            return points[i, 0] - points[i - 1, 0];
        }

        // - находит элемент матрицы С в главной диагонали.
        static public double Main_Diagonal(int i, double[,] points)
        {
            return (Get_Interval(i, points) + Get_Interval(i + 1, points)) / 3;
        }

        // - находит элемент матрицы С в верхней (нижней) диагонали.
        static public double Side_Diagonal(int i, double[,] points)
        {
            return Get_Interval(i, points) / 6;
        }

        // - находит вектор правой части.
        static public double Get_Vector(int i, double[,] points)
        {
            return ((points[i + 1, 1] - points[i, 1]) / Get_Interval(i + 1, points) - (points[i, 1] - points[i - 1, 1]) / Get_Interval(i, points));
        }

        // - решение СЛАУ методом прогонки.
        static public double[] Tridiagonal_Matrix_Algorithm(double[] lowerDiagonal, double[] b, double[] g, double[] d, int n)
        {
            double[] c = new double[n];
            double[] p = new double[n];
            double[] q = new double[n];
            for (int i = 0; i< n; i++)
            {
                if (i == 0)
                {
                    p[i] = -g[i] / b[i];
                    q[i] = d[i] / b[i];
                    continue;
                }
                p[i] = g[i] / (-b[i] - lowerDiagonal[i] * p[i - 1]);
                q[i] = (lowerDiagonal[i] * q[i - 1] - d[i]) / (-b[i] - lowerDiagonal[i] * p[i - 1]);
            }
            for (int i = n - 1; i>= 0; i--)
            {
                if (i == n - 1)
                {
                    c[i] = (lowerDiagonal[i] * q[i - 1] - d[i]) / (-b[i] - lowerDiagonal[i] * p[i - 1]);
                    continue;
                }
                c[i] = p[i] * c[i + 1] + q[i];
            }
            return c;
        }

        // - находит i-номер интервала, в котором лежит точка х.
        static public int Binary_Search(double[,] points, double x, int n)
        {
            int idx;
            if (x <= points[0, 0])
            {
                idx = 1;
            }
            else if (x >= points[n - 1, 0])
            {
                idx = n - 1;
            }
            else
            {
                int i = 0, j = n - 1;
                while (i + 1 < j)
                {
                    int k = i + (j - i) / 2;
                    if (x <= points[k, 0])
                    {
                        j = k;
                    }
                    else
                    {
                        i = k;
                    }
                }
                idx = j;
            }
            return idx;
        }

        // - суммирует элементы кубического сплайна.
        static public double Spline_Sum(int i, double[] Moments, double[,] points, double x)
        {
            double s1 = Moments[i - 1] * Math.Pow(points[i, 0] - x, 3) / (6 * Get_Interval(i, points));
            double s2 = Moments[i] * Math.Pow(x - points[i - 1, 0], 3) / (6 * Get_Interval(i, points));
            double s3 = (points[i - 1, 1] - Moments[i - 1] * Math.Pow(Get_Interval(i, points), 2) / 6) * (points[i, 0] - x) / Get_Interval(i, points);
            double s4 = (points[i, 1] - Moments[i] * Math.Pow(Get_Interval(i, points), 2) / 6) * (x - points[i - 1, 0]) / Get_Interval(i, points);
            return s1 + s2 + s3 + s4;
        }

        // - вычисляет кубический сплайн в точке x.
        static public double Cubic_Spline(double x, int n, double[,] points)
        {
            int i = 0;
            double s = 0;
            double[] nonZeroMoments = new double[n - 2];
            double[] Moments = new double[n];
            // нижняядиагональ
            double[] lowerDiagonal = new double[n - 2];
            // главнаядиагональ
            double[] mainDiagonal = new double[n - 2];
            // верхняядиагональ
            double[] upperDiagonal = new double[n - 2];
            // векторправыхчастей
            double[] d = new double[n - 2];
            lowerDiagonal[0] = upperDiagonal[n - 3] = 0;
            for (i = 1; i< n - 1; i++)
            {
                mainDiagonal[i - 1] = Main_Diagonal(i, points);
            }
            for (i = 2; i< n - 1; i++)
            {
                lowerDiagonal[i - 1] = Side_Diagonal(i, points);
            }
            for (i = 0; i< n - 3; i++)
            {
                upperDiagonal[i] = Side_Diagonal(i + 1, points);
            }
            for (i = 1; i< n - 2; i++)
            {
                d[i - 1] = Get_Vector(i, points);
            }
            nonZeroMoments = Tridiagonal_Matrix_Algorithm(lowerDiagonal, mainDiagonal, upperDiagonal, d, n - 2);
            for (i = 0; i< n - 2; i++)
            {
                Moments[i + 1] = nonZeroMoments[i];
            }
            Moments[0] = Moments[n - 1] = 0;
            i = Binary_Search(points, x, n);
            s = Spline_Sum(i, Moments, points, x);
            returns;
        }

        // - интегрирует функцию при помощи Элементарной формулы Симпсона.
        static public double Simpsons_Rule(double[,] points, int arr_size, double firstDerivative, int pow)
        {
            double x0 = points[0, 0];
            double y0 = points[0, 1];
            double x1 = points[arr_size, 0];
            double y1 = points[arr_size, 1];

            return (Math.Pow(y0, pow) + 4 * Math.Pow(Runge_Kutta_Method(x0, (x0 + x1) / 2, (x0 + x1) / 10, y0, firstDerivative), pow) + Math.Pow(y1, pow)) * (x1 - x0) / 6;
        }
    }
}
