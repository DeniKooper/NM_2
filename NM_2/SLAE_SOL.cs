using real = System.Double;
using System.IO;
using System;

namespace NM_2
{
    class SLAE_SOL
    {
        int n, m, end;
        real w, eps;
        real[] dia_0, dia_p1, dia_m1, dia_pm, dia_mm, dia_pm1, dia_mm1, dia_pm2, dia_mm2;
        real[] vector_b, vector_x0, vector_x1;
        real b_norm;
        real[] y;
        int l, k_b;


        public void get_sol_Jacobi(string path)
        {
            bool end_cycle = false;
            double nev;
            double x_rel = 0, x_norm = 0, mid = 0;
            for (int i = 0; i < n; i++)
                x_norm += (i + 1) * (i + 1);
            x_norm = Math.Sqrt(x_norm);
            FileStream output = new FileStream(path + "output_J.txt", FileMode.Create);
            StreamWriter os = new StreamWriter(output);
            int iter=0;


            for (int k = 1; k <= end && !end_cycle; k++)
            {

                double R_norm = 0;
                double temp_sum;

                base_sum();
                for (int i = 0; i < n; i++)
                {
                    R_norm += (y[i] - dia_0[i] * vector_x0[i]) * (y[i] - dia_0[i] * vector_x0[i]);
                }

                nev = Math.Sqrt(R_norm / b_norm);
                if (nev < eps) end_cycle = true;
                else
                {

                    for (int j = 0; j < n; j++)
                        vector_x1[j] = (1 - w) * vector_x0[j] + w * y[j] / dia_0[j];
                }


                real[] ch = vector_x0;
                vector_x0 = vector_x1;
                vector_x1 = ch;

                x_rel = 0;

                for (int i = 0; i < n; i++)
                    x_rel += (i + 1 - vector_x0[i]) * (i + 1 - vector_x0[i]);

                x_rel = Math.Sqrt(x_rel) / x_norm;
                mid += x_rel / nev;
                os.Write(k.ToString("5") + "\t" + (x_rel / nev).ToString("15") + Environment.NewLine);
                Console.Write("Iter: " + k + "\t Nev_sq: " + nev.ToString("4") + "\r");

                iter = k;
            }
            Console.Write(Environment.NewLine);
            os.Write("average: " + (mid / iter).ToString("15"));
        }

        public void get_sol_GaussSeidel(string path)
        {
            bool end_cycle = false;
            double nev;
            int end_d = n - m - 4;
            double x_rel = 0, x_norm = 0, mid = 0;
            for (int i = 0; i < n; i++)
                x_norm += (i + 1) * (i + 1);
            x_norm = Math.Sqrt(x_norm);
            FileStream output = new FileStream(path + "output_GS.txt", FileMode.Create);
            StreamWriter os = new StreamWriter(output);
            int iter = 0;

            for (int k = 1; k <= end && !end_cycle; k++)
            {

                double R_norm = 0;
                double temp_sum;

                base_sum();

                for (int i = 0; i < n; i++)
                {
                    R_norm += (y[i] - dia_0[i] * vector_x0[i]) * (y[i] - dia_0[i] * vector_x0[i]);
                }
                nev = Math.Sqrt(R_norm / b_norm);
                if (nev < eps) end_cycle = true;
                else
                {

                    for (int i = 0; i < end_d; i++)
                    {
                        vector_x1[i] = (1 - w) * vector_x0[i] + w * y[i] / dia_0[i];
                        y[i + 1] -= dia_m1[i] * (vector_x1[i] - vector_x0[i]);
                        y[i + m + 2] -= dia_mm[i] * (vector_x1[i] - vector_x0[i]);
                        y[i + m + 3] -= dia_mm1[i] * (vector_x1[i] - vector_x0[i]);
                        y[i + m + 4] -= dia_mm2[i] * (vector_x1[i] - vector_x0[i]);
                    }

                    vector_x1[end_d] = (1 - w) * vector_x0[end_d] + w * y[end_d] / dia_0[end_d];
                    y[end_d + 1] -= dia_m1[end_d] * (vector_x1[end_d] - vector_x0[end_d]);
                    y[end_d + m + 2] -= dia_mm[end_d] * (vector_x1[end_d] - vector_x0[end_d]);
                    y[end_d + m + 3] -= dia_mm1[end_d] * (vector_x1[end_d] - vector_x0[end_d]);

                    vector_x1[end_d + 1] = (1 - w) * vector_x0[end_d + 1] + w * y[end_d + 1] / dia_0[end_d + 1];
                    y[end_d + 1 + 1] -= dia_m1[end_d + 1] * (vector_x1[end_d + 1] - vector_x0[end_d + 1]);
                    y[end_d + 1 + m + 2] -= dia_mm[end_d + 1] * (vector_x1[end_d + 1] - vector_x0[end_d + 1]);

                    for (int i = end_d + 2; i < n; i++)
                        vector_x1[i] = (1 - w) * vector_x0[i] + w * y[i] / dia_0[i];
                    for (int i = end_d + 2; i < dia_m1.Length; i++)
                        y[i + 1] -= dia_m1[i] * (vector_x1[i] - vector_x0[i]);

                }


                real[] ch = vector_x0;
                vector_x0 = vector_x1;
                vector_x1 = ch;

                x_rel = 0;

                for (int i = 0; i < n; i++)
                    x_rel += (i + 1 - vector_x0[i]) * (i + 1 - vector_x0[i]);

                x_rel = Math.Sqrt(x_rel) / x_norm;
                mid += x_rel / nev;
                os.Write(k.ToString("5") + "\t" + (x_rel / nev).ToString("15") + Environment.NewLine);
                Console.Write("Iter: " + k + "\t Nev_sq: " + nev.ToString("4") + "\r");

                iter = k;
            }
            Console.Write(Environment.NewLine);
            os.Write("average: " + (mid / iter).ToString("15"));
        }


        void base_sum()
        {
            int[] ends = { n - m - 4, n - m - 3, n - m - 2, n - 1 };
            for (int i = 0; i < n; i++)
                y[i] = vector_b[i];

            for (int i = 0; i < ends[0]; i++)
                y[i] -= dia_pm2[i] * vector_x0[m + i + 4];
            for (int i = 0; i < ends[1]; i++)
                y[i] -= dia_pm1[i] * vector_x0[m + i + 3];
            for (int i = 0; i < ends[2]; i++)
                y[i] -= dia_pm[i] * vector_x0[m + i + 2];

            for (int i = 0; i < ends[3]; i++)
                y[i] -= dia_p1[i] * vector_x0[i + 1];
            for (int i = 0; i < ends[3]; i++)
                y[i + 1] -= dia_m1[i] * vector_x0[i];

            for (int i = 0; i < ends[2]; i++)
                y[m + 2 + i] -= dia_mm[i] * vector_x0[i];
            for (int i = 0; i < ends[1]; i++)
                y[m + 3 + i] -= dia_mm1[i] * vector_x0[i];
            for (int i = 0; i < ends[0]; i++)
                y[m + 4 + i] -= dia_mm2[i] * vector_x0[i];

        }


        public void blocks_r(string path)
        {
            calc_LU();
            bool end_cycle = false;
            int shift = 0;
            double nev = 0;
            double x_rel = 0, x_norm = 0, mid = 0;
            for (int i = 0; i < n; i++)
                x_norm += (i + 1) * (i + 1);
            x_norm = Math.Sqrt(x_norm);
            FileStream output = new FileStream(path + "output_B.txt", FileMode.Create);
            StreamWriter os = new StreamWriter(output);
            int iter = 0;


            for (int i = 0; i < end && !end_cycle; i++)
            {

                nev = nev_f() / Math.Sqrt(b_norm);

                if (nev < eps) end_cycle = true;
                if (!end_cycle)
                {
                    sol_sum_b();
                    for (int p = 0; p < n; p++)
                        y[p] *= w;
                    for (int j = 0; j < k_b; j++)
                    {
                        int[] ends = { (j + 1) * l, (j + 1) * l, (j + 1) * l };
                        if (ends[0] > n - m - 2) ends[0] = n - m - 2;
                        if (ends[1] > n - m - 3) ends[1] = n - m - 3;
                        if (ends[2] > n - m - 4) ends[2] = n - m - 4;

                        block_SLAE_sol(y, j);
                        shift = j * l;
                        for (int p = 0; p < l; p++, shift++)
                        {
                            vector_x1[shift] = (1 - w) * vector_x0[shift] + y[shift];
                        }
                        shift = (j + 1) * l;
                        if ((j + 1) * l < n)
                            y[shift - 1] += dia_m1[shift] * (vector_x0[shift] - vector_x1[shift]);

                        for (shift = j * l; shift < ends[0]; shift++)
                            y[shift + m + 2] += dia_mm[shift] * (vector_x0[shift] - vector_x1[shift]);

                        for (shift = j * l; shift < ends[1]; shift++)
                            y[shift + m + 3] += dia_mm[shift] * (vector_x0[shift] - vector_x1[shift]);

                        for (shift = j * l; shift < ends[2]; shift++)
                            y[shift + m + 4] += dia_mm[shift] * (vector_x0[shift] - vector_x1[shift]);
                    }

                    real[] ch = vector_x0;
                    vector_x0 = vector_x1;
                    vector_x1 = ch;

                    x_rel = 0;

                    for (int j = 0; j < n; j++)
                        x_rel += (j + 1 - vector_x0[j]) * (j + 1 - vector_x0[j]);

                    x_rel = Math.Sqrt(x_rel) / x_norm;
                    mid += x_rel / nev;
                    os.Write(i.ToString("5") + "\t" + (x_rel / nev).ToString("15") + Environment.NewLine);
                    Console.Write("Iter: " + i + "\t Nev_sq: " + nev.ToString("4") + "\r");

                    iter = i;
                }
                Console.Write(Environment.NewLine);
                os.Write("average: " + (mid / iter).ToString("15"));

            }
        }

        void block_SLAE_sol(real[] sol_vect, int b_n)
        {
            int shift = b_n * l;
            sol_vect[shift] /= dia_0[shift];
            shift++;
            for (int j = 1; j < l; j++, shift++)
                sol_vect[shift] = (sol_vect[shift] - dia_m1[shift - 1] * sol_vect[shift - 1]) / dia_0[shift];
            shift = l * (b_n + 1) - 2;
            for (int j = l - 2; j >= 0; j--, shift--)
                sol_vect[shift] -= dia_p1[shift] * sol_vect[shift + 1];
        }

        void calc_LU()
        {
            for (int i = 0; i < k_b; i++)
            {
                int shift = l * i;
                for (int j = 0; j < l - 1; j++, shift++)
                {
                    dia_p1[shift] /= dia_0[shift];
                    dia_0[shift + 1] -= dia_m1[shift] * dia_p1[shift];
                }
            }
        }

        void sol_sum_b()
        {
            int[] ends = { n - m - 4, n - m - 3, n - m - 2, n - 1 };
            for (int i = 0; i < n; i++)
                y[i] = vector_b[i];

            for (int i = 0; i < ends[0]; i++)
                y[i] -= dia_pm2[i] * vector_x0[m + i + 4];
            for (int i = 0; i < ends[1]; i++)
                y[i] -= dia_pm1[i] * vector_x0[m + i + 3];
            for (int i = 0; i < ends[2]; i++)
                y[i] -= dia_pm[i] * vector_x0[m + i + 2];

            for (int i = l - 1; i < ends[3]; i += l)
                y[i] -= dia_p1[i] * vector_x0[i + 1];
            for (int i = l - 1; i < ends[3]; i += l)
                y[i + 1] -= dia_m1[i] * vector_x0[i];

            for (int i = 0; i < ends[2]; i++)
                y[m + 2 + i] -= dia_mm[i] * vector_x0[i];
            for (int i = 0; i < ends[1]; i++)
                y[m + 3 + i] -= dia_mm1[i] * vector_x0[i];
            for (int i = 0; i < ends[0]; i++)
                y[m + 4 + i] -= dia_mm2[i] * vector_x0[i];
        }

        public SLAE_SOL(StreamReader m_file, StreamReader v_file, StreamReader b_file){
            enter_data(m_file, v_file, b_file);
        }

        void enter_data(StreamReader m_file, StreamReader v_file, StreamReader b_file)
        {
            String[] temp = m_file.ReadLine().Split(' ');

            n = Convert.ToInt32(temp[0]);
            m = Convert.ToInt32(temp[1]);
            w = Convert.ToInt32(temp[2]);
            eps = Convert.ToDouble(temp[3]);
            end = Convert.ToInt32(temp[4]);
            l = Convert.ToInt32(temp[5]);
            
            k_b = n / l;

            dia_0 = new double[n];
            dia_p1 = new double[n - 1];
            dia_m1 = new double[n - 1];
            dia_pm = new double[n - m - 2];
            dia_mm = new double[n - m - 2];
            dia_pm1 = new double[n - m - 3];
            dia_mm1 = new double[n - m - 3];
            dia_pm2 = new double[n - m - 4];
            dia_mm2 = new double[n - m - 4];

            vector_b = new double[n];
            vector_x0 = new double[n];
            vector_x1 = new double[n];

            y = new double[n];

            readVector(m_file, n-m-4, dia_pm2);
            readVector(m_file, n-m-3, dia_pm1);
            readVector(m_file, n-m-2, dia_pm);

            readVector(m_file, n-1, dia_p1);
            readVector(m_file, n, dia_0);
            readVector(m_file, n-1, dia_m1);

            readVector(m_file, n-m-2, dia_mm);
            readVector(m_file, n-m-3, dia_mm1);
            readVector(m_file, n-m-4, dia_mm2);

            b_norm = 0;
            temp = v_file.ReadLine().Split(' ');
            for (int i = 0; i < n; i++)
            {
                vector_b[i] = Convert.ToDouble(temp[i]);
                b_norm += vector_b[i] * vector_b[i];
            }

            readVector(b_file, n, vector_x0);            
        }

        void readVector(StreamReader fs, int n, real[] v)
        {
            String[] temp = fs.ReadLine().Split(' ');
            for (int i = 0; i < n; i++)
                v[i] = Convert.ToDouble(temp[i]);
        }

        public void output_solution(StreamWriter out_f)
        {
            for (int i = 0; i < n; i++)
                out_f.Write(vector_x0[i].ToString("E15") + Environment.NewLine);
        }

        double nev_f()
        {
            sol_sum_b();
            double nev_norm = 0;
            double s;
            int shift;

            for (int j = 0; j < k_b; j++)
            {
                shift = j * l;
                s = y[shift] - dia_0[shift] * vector_x0[shift] - dia_p1[shift] * dia_0[shift] * vector_x0[shift + 1];
                nev_norm += s * s;
                shift++;
                for (int i = 1; i < l - 1; i++, shift++)
                {
                    s = y[shift] - (dia_0[shift] + dia_m1[shift - 1] * dia_p1[shift - 1]) * vector_x0[shift] - dia_m1[shift - 1] * vector_x0[shift - 1] - dia_p1[shift] * dia_0[shift] * vector_x0[shift + 1];
                    nev_norm += s * s;
                }
                s = y[shift] - (dia_0[shift] + dia_m1[shift - 1] * dia_p1[shift - 1]) * vector_x0[shift] - dia_m1[shift - 1] * vector_x0[shift - 1];
                nev_norm += s * s;

            }
            nev_norm = Math.Sqrt(nev_norm);
            return nev_norm;
        }

    }
}


 