package ppa_program;
import java.io.*;

/** Readers for PWM matrices, pentamer probabilities and NUE hexamer weights. */
public class DataLoader {

    //*********************************************************************
    // nue weight initilaize
    //*******************************************************************
    public static void nueWeightInitialize(double[][] hexa, String file) {
        try {
            BufferedReader b = new BufferedReader(new FileReader(file));
          //  BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Users/Nadeem/Desktop/output.txt"));

            String line = new String();
            String num = new String();
            String use = new String();
            String dse = new String();
            int i = 0, index = 0, pentamerNum = 0;
            char alpha;


            line = b.readLine();
            while (line != null) {
                i++;
                index = 0;
                num = new String();
                use = new String();
                dse = new String();
                char[] arr = new char[line.length()];
                line.getChars(0, (line.length()), arr, 0);
                while (index < 6) {
                    alpha = arr[index++];
                    num = num + MotifScoring.getNum(alpha);
                }
                index += 5;
                while (line.charAt(index) != ' ') {
                    use = use + arr[index++];
                }
                while (line.charAt(index) == ' ')
                    index++;
                while (index != line.length()) {
                    dse = dse + arr[index++];
                }
                pentamerNum = Integer.parseInt(num);

                hexa[pentamerNum][0] = Double.parseDouble(use);
                hexa[pentamerNum][1] = Double.parseDouble(dse);
              //  System.out.println("" + pentamerNum + " Count = " + hexa[pentamerNum][0] + " Weiht = " + hexa[pentamerNum][1]);
               // bw.write("" + pentamerNum + " Count = " + hexa[pentamerNum][0] + " Weiht = " + hexa[pentamerNum][1]);
               // bw.newLine();
                line = b.readLine();
            }
          //  bw.close();

        } catch (IOException e) {
            System.out.println("Exception : " + e.getMessage());
        }

    }    // end of nueWeightInitilize


    //*********************************************************************
    // penta Initilization
    //*******************************************************************
    public static void pentaInitialize(double[][] hexa, String file) {
        try {
            BufferedReader b = new BufferedReader(new FileReader(file));
          //  BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Users/Nadeem/Desktop/output.txt"));

            String line = new String();
            String num = new String();
            String nue = new String();
            String fue = new String();
            String dse = new String();
            int i = 0, index = 0, pentamerNum = 0;
            char alpha;


            line = b.readLine();
            while (line != null) {
                i++;
                index = 0;
                num = new String();
                nue = new String();
                fue = new String();
                dse = new String();
                char[] arr = new char[line.length()];
                line.getChars(0, (line.length()), arr, 0);
                while (index < 5) {
                    alpha = arr[index++];
                    num = num + MotifScoring.getNum(alpha);
                }
                index += 3;
                while (line.charAt(index) != ' ') {
                    dse = dse + arr[index++];
                }
                while (line.charAt(index) == ' ')
                    index++;
                while (line.charAt(index) != ' ') {
                    nue = nue + arr[index++];
                }
                while (line.charAt(index) == ' ')
                    index++;
                while (index != line.length()) {
                    fue = fue + arr[index++];
                }
                pentamerNum = Integer.parseInt(num);

                hexa[pentamerNum][0] = Double.parseDouble(dse);
                hexa[pentamerNum][1] = Double.parseDouble(nue);
                hexa[pentamerNum][2] = Double.parseDouble(fue);
              //  System.out.println("" + pentamerNum + " Count = " + hexa[pentamerNum][0] + " Weiht = " + hexa[pentamerNum][1]);
             //   bw.write("" + pentamerNum + "   " + hexa[pentamerNum][0] + "     " + hexa[pentamerNum][1] + "     " + hexa[pentamerNum][2]);
            //    bw.newLine();
                line = b.readLine();
            }
          //  bw.close();

        } catch (IOException e) {
            System.out.println("Exception : " + e.getMessage());
        }

    }    // end of pentaInitialize


    //*********************************************************************
    // Tri Initilization
    //*******************************************************************
    public static void triInitializing(double[] tri, File file) {
        try {
            BufferedReader b = new BufferedReader(new FileReader(file));

            String line = new String();
            String num = new String();
            String use = new String();
            String dse = new String();
            int i = 0, index = 0, pentamerNum = 0;
            char alpha;


            line = b.readLine();
            while (line != null) {
                i++;
                index = 0;
                num = new String();
                use = new String();
                dse = new String();
                char[] arr = new char[line.length()];
                line.getChars(0, (line.length()), arr, 0);
                while (index < 3) {
                    alpha = arr[index++];
                    num = num + MotifScoring.getNum(alpha);
                }
                index += 3;

                while (index != line.length()) {
                    dse = dse + arr[index++];
                }
                pentamerNum = Integer.parseInt(num);

                tri[pentamerNum] = Double.parseDouble(dse);
                line = b.readLine();
            }

        } catch (IOException e) {
            System.out.println("Exception : " + e.getMessage());
        }

    }    // end of tri initialization


    //*********************************************************
    public static void readMatrixFromFile(double[][] matrix,String filePath){

    	// Initialization
		String line = "",score = "";
		int countLine = 0;

		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
			line = br.readLine();
			line = br.readLine();
			line = br.readLine();
			while(line != null){
				int index = 0, ele = 0;

				countLine++;
				while(line.charAt(index) == ' ')
					index++;
				while(index != line.length()){


					if(line.charAt(index) == ' '){
						matrix[countLine - 1][ele++] = Double.parseDouble(score);
					//	System.out.println(" " + countLine + " : " + matrix[countLine - 1][ele - 1]);
						score = "";
						while(line.charAt(index) == ' ')
							index++;
					}
					else{
						score += line.charAt(index);
						index++;
					}
				}
				if(index == line.length()){
					matrix[countLine - 1][ele++] = Double.parseDouble(score);score = "";}
				line = br.readLine();
			} // end of main while
                        br.close();
		}  // end of try block

		catch(IOException e){
			System.out.println(" Exception : " + e.getMessage() );
		}

    }


}
