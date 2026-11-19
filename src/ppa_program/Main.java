package ppa_program;
import java.io.*;

/**
 * PPARF: plant poly(A) site prediction.
 *
 * Usage: java -cp out ppa_program.Main [matrixDir] [input.fasta] [output.txt]
 * Defaults (original setup): C:\ppa_matrices  D:\Loke_8k_train.seq  D:\pparf_output.txt
 */
public class Main {

    public static void main(String[] args) {
       // Command-line paths (fall back to the original hardcoded locations)
		   String matrixDir = args.length > 0 ? args[0] : "C:\\ppa_matrices";
		   String inputPath = args.length > 1 ? args[1] : "D:\\Loke_8k_train.seq";
		   String outputPath = args.length > 2 ? args[2] : "D:\\pparf_output.txt";

       // Variables & Initializations
		   String line = "",idLine = "",seq = "";
		   int analyzed = 0, count = 0;

		  // String csPath = "/home/Nadeem/Desktop/ppa_EM_matrix/cs.txt";
		   String nuePath = new File(matrixDir, "nue_ji_final.txt").getPath();
		   String cePath = new File(matrixDir, "ce_m2.res").getPath();
		   String fueShortPath = new File(matrixDir, "fue_short.res").getPath();
		   String fueLongPath = new File(matrixDir, "fue_long.res").getPath();
		   String ndePath = new File(matrixDir, "nde.res").getPath(); // in matrices folder
		   //String pentaFilePath = "C:/Users/Nadeem/Desktop/ppa/pentamer_prob_1stapproach(mix).txt";
		   String pentaFilePath = new File(matrixDir, "pentamer_prob_cds.txt").getPath(); // oligomat - >50
		   String nueWeightFilePath = new File(matrixDir, "weighted_hexamers_NUE_new.txt").getPath(); // oligomat -> 60


		// Arrays
		   double[][] signal = new double[100000][2]; // 0->position, 1-> score
		   double[][] nueMatrix = new double[4][6];
		   double[][] ceMatrix = new double[4][8];
		   double[][] csMatrix = {{0.2657,0.2294,0.2053,0.2003,0.2019,0.2217,0.2327,0.236,0.1533,0.0589,0.733,0.4697,0.3053,0.2771,0.287,0.2499,0.2643,0.2474,0.2253,0.2171},{0.1383,0.1361,0.1216,0.1221,0.1391,0.1486,0.1503,0.1954,0.2161,0.3153,0.0606,0.1074,0.1359,0.1606,0.1649,0.1721,0.159,0.1563,0.1531,0.1671},{0.1584,0.1509,0.1403,0.1249,0.1353,0.1597,0.1663,0.1607,0.1423,0.139,0.0433,0.0817,0.1233,0.1539,0.1531,0.1636,0.1714,0.1869,0.1883,0.1813},{0.4376,0.4836,0.5329,0.5527,0.5237,0.47,0.4507,0.4079,0.4883,0.4869,0.1631,0.3411,0.4356,0.4084,0.395,0.4144,0.4053,0.4094,0.4333,0.4344}};
		   double[] ncDis = {0.1349,0.2857,0.3324,0.2470,0,0,0};  // -30:-5
		   double[] fcDis = {0.1086,0.1056,0.1130,0.1081,0.1186,0.1090,0.094,0.093,0.0867,0.0104,0};
		   double[][] fueShortMatrix = new double[4][6];
		   double[][] fueLongMatrix = new double[4][10];
		   double[][] ndeMatrix = new double[4][10];
		   double[][] hexa = new double[44445][3]; // 0 => pentaNDE, 1 => pentaNUE, 2 => pentaFUE
		   double[][] nueWeight = new double[444445][2]; // 0 => counts, 1 => weights

		   double weight = 0;
		   String[] nueResult = new String[3];     // 0 => motif, 1 => score, 2 => indexInseq
		   String[] ceResult = new String[3];       // 0 => motif, 1 => score, 2 => indexInseq
		   String[] fueShortResult = new String[3];      // 0 => motif, 1 => score, 2 => indexInseq
		   String[] fueLongResult = new String[3];
		   String[] ndeResult = new String[3];
		   //String[] csResult = new String[3];
		   double[] total = new double[13];
		   int arrayIndex = 0,totalPos= 0;
		   for(int i = 0; i < 13; i++)
			   total[i] = 0;
		   for (int i = 0; i < 44445; i++)
		   {
			   hexa[i][0] = 0; hexa[i][1] = 0; hexa[i][2] = 0;
		   }

		   DataLoader.readMatrixFromFile(ceMatrix, cePath);
		   DataLoader.readMatrixFromFile(fueShortMatrix, fueShortPath);
		   DataLoader.readMatrixFromFile(fueLongMatrix, fueLongPath);
		   DataLoader.readMatrixFromFile(ndeMatrix, ndePath);
		   DataLoader.readMatrixFromFile(nueMatrix, nuePath);
		   DataLoader.nueWeightInitialize(nueWeight, nueWeightFilePath);
		   DataLoader.pentaInitialize(hexa, pentaFilePath);

		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));//Leave_one_out\\results_leave_weight.txt"));
			//BufferedWriter bw1 = new BufferedWriter(new FileWriter("C:\\Leave_one_out\\results1.txt"));
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(inputPath)));
		    String w = "NUE     FUElong     NuePentaScore     DownPentaScore     FuePentaScore     CS     NUECSdis";//     FUEratio     NUEratio     CEratio";
		    bw.write(w);
		    bw.newLine();
		    line = br.readLine();
		    while(line != null){
		    	if(line.startsWith(">")){
		    		bw.write(line);
		    		bw.newLine();
		    		analyzed++;
		    		seq = "";
		    		idLine = line;
		    		line = br.readLine();
		    		while(!(line.startsWith(">"))){
		    			seq += line;
		    			line = br.readLine();
		    			if(line == null)
		    				break;
		    		}   // end of get seq while
		    		seq = seq.toUpperCase();
					int loop = 150, status = 0,pos = 0;
					arrayIndex = 0;
					while (loop < seq.length() - 61)
					{
						String hexaNUE = "", hexaFUE = "", hexaDOWN = "", nueRegion = "", fueShortRegion = "", fueLongRegion = "", ceRegion = "", ndeRegion = "", csRegion = "", fueRatioRegion = "", nueRatioRegion = "", ceRatioRegion = "";

						for (int i = (loop - 50); i <= (loop - 1); i++)
							hexaNUE += seq.charAt(i);
						for (int i = loop; i <= (loop + 59); i++)
							hexaDOWN += seq.charAt(i);
						for (int i = (loop - 150); i <= (loop - 150+99); i++)
							hexaFUE += seq.charAt(i);
					//	for (int i = (loop - 40); i <= (loop - 1); i++)
						//	nueRegion += seq.charAt(i);
						for (int i = (loop - 100); i <= (loop - 41); i++)
							fueShortRegion += seq.charAt(i);
						for (int i = (loop - 140); i <= (loop - 41); i++)
							fueLongRegion += seq.charAt(i);
						for (int i = (loop - 20); i <= (loop + 19); i++)
							ceRegion += seq.charAt(i);
						for (int i = (loop - 10); i <= (loop + 9); i++)
							csRegion += seq.charAt(i);
						for (int i = (loop - 41); i <= (loop - 11); i++)
							nueRatioRegion += seq.charAt(i);
						for (int i = (loop - 150); i <= (loop - 36); i++)
							fueRatioRegion += seq.charAt(i);
						for (int i = (loop - 10); i <= (loop + 19); i++)
							ceRatioRegion += seq.charAt(i);
						for(int i = (loop - 30); i <= (loop - 6); i++)
							nueRegion += seq.charAt(i);



						MotifScoring.NUE(nueRegion, nueMatrix, nueResult, nuePath, 6);
						// weight portion
						String nueNum = "", motif = nueResult[0];
						nueNum = nueNum + MotifScoring.getNum(motif.charAt(0)) + MotifScoring.getNum(motif.charAt(1)) + MotifScoring.getNum(motif.charAt(2)) + MotifScoring.getNum(motif.charAt(3)) + MotifScoring.getNum(motif.charAt(4)) + MotifScoring.getNum(motif.charAt(5));
						weight = nueWeight[Integer.parseInt(nueNum)][1];

						MotifScoring.FUE(fueLongRegion, fueLongMatrix, fueLongResult, fueLongPath, 10);

		                double fueScore =(Double.parseDouble(fueLongResult[1]));// * MotifScoring.disFUE(fcDis,Integer.parseInt(fueLongResult[2]),Double.parseDouble(fueLongResult[1]));
						total[2] += Double.parseDouble(fueLongResult[1]);

						double csScore = MotifScoring.csScore(csMatrix, csRegion);

					//	if((Double.parseDouble(nueResult[1]) >= 1.912) && (fueScore >=3.086) && (csScore >= 4.1957)){
							double score = ((Double.parseDouble(nueResult[1])));
							total[0] += score;

							double nueHexaScore = MotifScoring.hexaScore(hexaNUE, hexa, 1);
							double downHexaScore = MotifScoring.hexaScore(hexaDOWN, hexa, 0);
							double fueHexaScore = MotifScoring.hexaScore(hexaFUE, hexa, 2);


							double distance = MotifScoring.nuecsDis(ncDis, nueResult);
							double fueRatio = MotifScoring.ratioUA(fueRatioRegion);
							double ceRatio = MotifScoring.ratioUA(ceRatioRegion);
							double nueRatio = MotifScoring.ratioAU(nueRatioRegion);


							// cds coefficients
							double integral = ((score/6)*4.515) + ((fueScore/10) * 2.67) + ((nueHexaScore / 46)*28.019) + ((downHexaScore / 56)*28.689) + ((fueHexaScore / 96)*45.846) + ((csScore / 20)*39.541) + ((distance)*3.816)+(weight*0.267);

							if(integral >= 68.043){
								pos = 1;
								if(((loop+1)>=140) && ((loop+1)<=160))
									status = 1;
								signal[arrayIndex][0] = loop+1;
								signal[arrayIndex++][1] = integral;
							}
						//} // end of min score IF

						loop++;
					}

		            for(int k = 0; k < arrayIndex; k++){
		            	bw.write("Pos: " + signal[k][0] + "     Score: " + signal[k][1]);
		            	bw.newLine();
		            }
		            if(status == 1)
		            	count++;
		            if(pos == 1)
		            	totalPos++;
		            pos = 0;
		            status = 0;

		    	} // end of ">" if
		    	else
		    		line = br.readLine();
		    } // End of while

			System.out.println(" Analyzed : " + analyzed + " NUE = " + MotifScoring.roundTwoFour("" + (total[0] / analyzed)) + " FUEshort = " + MotifScoring.roundTwoFour("" + (total[1] / analyzed)) + " FUElong = " + MotifScoring.roundTwoFour("" + (total[2] / analyzed)) + " CE = " + MotifScoring.roundTwoFour("" + (total[3] / analyzed)) + " NDE = " + MotifScoring.roundTwoFour("" + (total[4] / analyzed)) + " NUEHexa = " + MotifScoring.roundTwoFour("" + (total[5] / analyzed)) + " DownHexa = " + MotifScoring.roundTwoFour("" + (total[6] / analyzed)) + " FUEHexa = " + MotifScoring.roundTwoFour("" + (total[7] / analyzed)) + " CS = " + MotifScoring.roundTwoFour("" + (total[8] / analyzed)) + " distance = " + MotifScoring.roundTwoFour("" + (total[9] / analyzed)) + " fueRatio = " + MotifScoring.roundTwoFour("" + (total[10] / analyzed)) + " ceRatio = " + MotifScoring.roundTwoFour("" + (total[11] / analyzed)) + " nueRatio = " + MotifScoring.roundTwoFour("" + (total[12] / analyzed)));
		    //bw1.close();
		    br.close();
		    bw.close();
		    System.out.println("count : " + count + " totalPos = " + totalPos);
		} // end of try block
		catch(IOException e){
			System.out.println("Exception : " + e.getMessage());
		}
		finally{
			System.out.println("Processs Completed : ");
		}

	} // end of main method
}
