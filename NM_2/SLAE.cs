using System;
using System.IO;
using real = System.Double;

namespace NM_2
{
    class SLAE
    {
        real[] d0;
        real[] du1;
        real[] du2;
        real[] du3;
        real[] du4;
        real[] dl1;
        real[] dl2;
        real[] dl3;
        real[] dl4;
        public real[] x0;//начальное приближение
        public real[] x1;//текущий вектор
        public real[] b; //вектор правой части
        int[] I; //вектор смещений
                 //размеры диагоналей матрицы
        public int[] len;
        int sizeBlock;
        //точность решения
        real eps;
        public double w;
        //макисмальное количество интераций
        int maxiter;
        int metod; //0-Якоби, 1-Гаусс Зейдель

        String path;
        public SLAE(String path)
        {
            this.path = path;
            Input();
        }

        public void Input()
        {
            Console.WriteLine("Preparing to read a input.txt");
            FileStream input = new FileStream(path + "input.txt", FileMode.Open);
            StreamReader inputStream = new StreamReader(input);
            Console.WriteLine("Preparing to read a inf.txt");
            FileStream inf = new FileStream(path + "inf.txt", FileMode.Open);
            StreamReader infStream = new StreamReader(inf);
            Console.WriteLine("Preparing to read a vector.txt");
            FileStream vector = new FileStream(path + "vector.txt", FileMode.Open);
            StreamReader vectorStream = new StreamReader(vector);

            len = new int[5];

            readIntVector(inputStream, len);

            I = new int[5];
            d0 = new real[len[0]];
            du1 = new real[len[1]];
            du2 = new real[len[2]];
            du3 = new real[len[3]];
            du4 = new real[len[4]];
            dl1 = new real[len[1]];
            dl2 = new real[len[2]];
            dl3 = new real[len[3]];
            dl4 = new real[len[4]];
            x0 = new real[len[0]];
            x1 = new real[len[0]];
            b = new real[len[0]];

            readIntVector(infStream, I);
            eps = Convert.ToDouble(infStream.ReadLine());
            w = Convert.ToDouble(infStream.ReadLine());
            maxiter = Convert.ToInt32(infStream.ReadLine());

            metod = Convert.ToInt32(infStream.ReadLine());
            sizeBlock = Convert.ToInt32(infStream.ReadLine());

            if (metod == 0)
                readVector(infStream, len[0], x0);
            else
                readVector(infStream, len[0], x1);

            readVector(inputStream, len[0], d0);
            readVector(inputStream, len[1], dl1);
            readVector(inputStream, len[2], dl2);
            readVector(inputStream, len[3], dl3);
            readVector(inputStream, len[4], dl4);
            readVector(inputStream, len[1], du1);
            readVector(inputStream, len[2], du2);
            readVector(inputStream, len[3], du3);
            readVector(inputStream, len[4], du4);
            readVector(vectorStream, len[0], b);

            inputStream.Close();
            infStream.Close();
            vectorStream.Close();
        }
        
        void readIntVector(StreamReader fs, int[] v)
        {
            String[] temp = fs.ReadLine().Split(' ');
            for (int i = 0; i < 5; i++)
                v[i] = Convert.ToInt32(temp[i]);
        }

        void readVector(StreamReader fs, int n, real[] v)
        {
            String[] temp = fs.ReadLine().Split(' ');
            for (int i = 0; i < n; i++)
                v[i] = Convert.ToDouble(temp[i]);
        }

        //умножение строки на вектор
        double mult(int i, real[] x)
        {
            double sum = 0;
            sum = d0[i] * x[i];
            //обрабатываем нижние диагонали
            if (i >= I[1])
                sum += dl1[i - I[1]] * x[i - I[1]];
            if (i >= I[2])
                sum += dl2[i - I[2]] * x[i - I[2]];
            if (i >= I[3])
                sum += dl3[i - I[3]] * x[i - I[3]];
            if (i >= I[4])
                sum += dl4[i - I[4]] * x[i - I[4]];
            //обабатываем верхние диагонали
            if (i < len[1])
                sum += du1[i] * x[i + I[1]];
            if (i < len[2])
                sum += du2[i] * x[i + I[2]];
            if (i < len[3])
                sum += du3[i] * x[i + I[3]];
            if (i < len[4])
                sum += du4[i] * x[i + I[4]];
            return sum;
        }

        double multAx(int i, real[] x, int j1, int j2)
        {
            double sum = 0;
            //sum = d0[i]*x[i];
            //обрабатываем нижние диагонали
            if (i >= I[1] && ((i <= j1 || i > j2) && i != len[0] - 1))
                sum += dl1[i - I[1]] * x[i - I[1]];
            if (i >= I[2])
                sum += dl2[i - I[2]] * x[i - I[2]];
            if (i >= I[3])
                sum += dl3[i - I[3]] * x[i - I[3]];
            if (i >= I[4])
                sum += dl4[i - I[4]] * x[i - I[4]];
            //обабатываем верхние диагонали
            if (i < len[1] && ((i < j1 || i >= j2) && i != 0))
                sum += du1[i] * x[i + I[1]];
            if (i < len[2])
                sum += du2[i] * x[i + I[2]];
            if (i < len[3])
                sum += du3[i] * x[i + I[3]];
            if (i < len[4])
                sum += du4[i] * x[i + I[4]];
            return sum;
        }

        double cond_nev()
        {
            double cond = 0;
            //погрешность 
            double p_x = 0;
            double norm_x = 0;
            for (int i = 0; i < len[0]; i++)
            {
                p_x += Math.Pow(x1[i] - (i + 1), 2);
                norm_x += Math.Pow(i + 1, 2);
            }
            norm_x = Math.Sqrt(p_x / norm_x);
            //невязка 
            double Ax = 0;
            double p_F = 0;
            double norm_F = 0;
            for (int i = 0; i < len[0]; i++)
            {
                Ax = mult(i, x1);
                p_F += Math.Pow(b[i] - Ax, 2);
                norm_F += Math.Pow(b[i], 2);
            }
            norm_F = Math.Sqrt(p_F / norm_F);
            cond = norm_x / norm_F;
            return cond;
        }

        public double condA(real[] x, double normB)
        {
            double cond = 0;
            double xNorm = 0;
            double sum = 0;
            double pog = 0;
            for (int i = 0; i < len[0]; i++)
            {
                sum += Math.Pow(x[i] - (i + 1), 2);
                xNorm += Math.Pow(i + 1, 2);
            }
            pog = Math.Sqrt(sum / xNorm);

            //невязка 
            sum = 0;
            double sumNev = 0;
            double nev2 = 0;
            double nev = 0;
            for (int i = 0; i < len[0]; i++)
            {
                sum = mult(i, x1);
                sumNev = b[i] - sum;
                nev2 += Math.Pow(sumNev, 2); //невзяка 
            }
            nev = Math.Sqrt(nev2 / normB);   //относительная невязка
            cond = pog / nev;
            return cond;
        }

        public void clearV(real[] x)
        {
            for (int i = 0; i < len[0]; i++)
                x[i] = 0;
        }

        public void makeLU()
        {
            int colBlock = len[0] / sizeBlock; //количество блоков
            int offset = 0;                 //смещение
            for (int i = 0; i < colBlock; i++)
            {
                offset = sizeBlock * i;
                for (int j = 0; j < sizeBlock - 1; j++, offset++)
                {
                    du1[offset] /= d0[offset];              //обработка верхнего треугольника
                    d0[offset + 1] -= du1[offset] * dl1[offset];//обработка главной диагонали 
                }
            }
        }

        void multLY(int i, real[] f)
        {
            //x0 = y
            f[i] = f[i] / d0[i];
            int oldi = i;
            for (i = i + 1; i < sizeBlock + oldi; i++)
            {
                double sum = 0;
                sum = f[i - 1] * dl1[i - 1];
                f[i] = (f[i] - sum) / d0[i];
            }
        }

        void multUX(int i, real[] f)
        {
            f[i] = f[i];
            int oldi = i;
            for (i = i - 1; i >= oldi - sizeBlock + 1; i--)
            {
                f[i] = f[i] - f[i + 1] * du1[i];
            }
        }

        public void blockRelax(StreamWriter os)
        {
            /*std::ofstream os("output.txt");*/
            int outB = 0, inB = 0, j1 = 0, j2 = 0;
            double sum = 0;
            int iter = 0;
            double nev = 0, nevSum = 0, nevSum2 = 0;
            bool exit = false;
            real[] R;
            R = new real[len[0]];
            real[] Y;
            Y = new real[len[0]];
            int colBlock = len[0] / sizeBlock; //количество блоков
            for (int k = 0; k < maxiter && !exit; k++)
            {

                //обрабатываем по блокам
                for (int i = 0; i < len[0]; i += sizeBlock)
                {
                    j1 = i;
                    j2 = i + sizeBlock - 1;
                    //внутри блока
                    for (int ib = i; ib < sizeBlock + i; ib++)//по строкам
                    {
                        sum = multAx(ib, x1, j1, j2);
                        R[ib] = w * (b[ib] - sum);
                    }
                    multLY(j1, R);
                    multUX(j2, R);
                    //x1 = R;
                    for (int r = i; r < sizeBlock + i; r++)
                        x1[r] = R[r] + (1 - w) * x0[r];
                }

                iter = k;
                nevSum = 0;
                nevSum2 = 0;
                for (int v = 0; v < len[0]; v++)
                {
                    nevSum += Math.Pow((x1[v] - x0[v]), 2);
                    nevSum2 += Math.Pow(x1[v], 2);
                }
                nev = Math.Sqrt(nevSum / nevSum2);
                if (nev < eps)
                {
                    exit = true;
                    iter = k;
                    Console.WriteLine("Выход по невязке\n");
                    break;
                }

                x0 = x1;
            }
            nevSum = 0;
            nevSum2 = 0;
            for (int v = 0; v < len[0]; v++)
            {
                nevSum += Math.Pow((x1[v] - x0[v]), 2);
                nevSum2 += Math.Pow(x0[v], 2);
            }
            nev = Math.Sqrt(nevSum / nevSum2);
            os.Write("-----------------------------" + Environment.NewLine);
            os.Write(w.ToString("3") + "\t"+ iter + Environment.NewLine);
            for (int i = 0; i < len[0]; i++)
                os.Write(x1[i].ToString("18")+ Environment.NewLine);
        }
    }
}
