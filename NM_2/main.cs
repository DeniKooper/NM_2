using System;
using System.IO;

namespace NM_2
{
    class main
    {
        static void Main(string[] args)
        {
            SLAE_SOL our_S;
            string path = "C:\\Users\\Kuhti\\source\\repos\\NM_2\\NM_2\\test\\";

            FileStream f1 = new FileStream(path + "b_mp.txt", FileMode.Open);
            StreamReader matrix = new StreamReader(f1);
            FileStream f2 = new FileStream(path + "b_vp.txt", FileMode.Open);
            StreamReader v_b = new StreamReader(f2);
            FileStream f3 = new FileStream(path + "x0.txt", FileMode.Open);
            StreamReader v_x0 = new StreamReader(f3);

            our_S= new SLAE_SOL(matrix, v_b, v_x0);

            f1.Close();
            f2.Close();
            f3.Close();

            FileStream f4 = new FileStream(path + "output.txt", FileMode.Create);
            StreamWriter output = new StreamWriter(f4);
            //our_S.get_sol_Jacobi(path);
            //our_S.get_sol_GaussSeidel(path);
            //our_S.blocks_r(path);
            our_S.output_solution(output);
            f4.Close();
        }
    }
}
