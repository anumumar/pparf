package ppa_program;

/** PWM motif search, pentamer scoring, CS profile and distance features. */
public class MotifScoring {

// ***********************************************
	public static void NUE(String seq, double[][] scoreMatrix , String[] tataResult, String filePath, int len){


        int start = 0, end = seq.length() - 1, bestIndex = 0;
		len = len - 1;
        double score = 0,bestScore = 0;
        String motif = "",newMotif = "";
        while(start <= end){
       	 score = 0;
       	 newMotif = "";
       	 for(int i = start; i <= (start + len); i++){
       		 if(i < seq.length()){
       			 newMotif+=seq.charAt(i);
       			if(seq.charAt(i) == 'A')
       				score += scoreMatrix[0][i - start];
       			else if(seq.charAt(i) == 'C')
       				score += scoreMatrix[1][i - start];
       			else if(seq.charAt(i) == 'G')
       				score += scoreMatrix[2][i - start];
       			else{
       				if(seq.charAt(i) == 'T')
       			     	score += scoreMatrix[3][i - start];
       			}
       		 } // end of IF for sequence length
       	 } // end of for


       	 if(bestScore == 0){
       		 bestScore = score;
       		 motif = newMotif;
       		 bestIndex = start + 1;
       	 }
       	 else{
       		 if((score > bestScore) && (newMotif.length() == 6)){
       			 bestScore = score;
       			 motif = newMotif;
       			 bestIndex = start + 1;
       		 }
       	 }
       	 start++;
        } // end of while
	    tataResult[0] = motif;
	    tataResult[1] = "" + bestScore;
	    tataResult[2] = "" + bestIndex;
	}  // end of TATA Method

	// ***********************************************
	public static void FUE(String seq, double[][] scoreMatrix , String[] ypResult, String filePath, int len){


		//double upHexaScore = hexaScore(hexaUP, hexa, 0);
        int start = 0, end = seq.length() - 1, bestIndex = 0;
        len = len - 1;
        double score = 0,bestScore = 0;
        String motif = "",newMotif = "";
        while(start <= end){
       	 score = 0;
       	 newMotif = "";
       	 for(int i = start; i <= (start + len); i++){
       		 if(i < seq.length()){
       			 newMotif+=seq.charAt(i);
       			if(seq.charAt(i) == 'A')
       				score += scoreMatrix[0][i - start];
       			else if(seq.charAt(i) == 'C')
       				score += scoreMatrix[1][i - start];
       			else if(seq.charAt(i) == 'G')
       				score += scoreMatrix[2][i - start];
       			else{
       				if(seq.charAt(i) == 'T')
       			     	score += scoreMatrix[3][i - start];
       			}
       		 } // end of IF for sequence length
       	 } // end of for


       	 if(bestScore == 0){
       		 bestScore = score;
       		 motif = newMotif;
       		 bestIndex = start + 1;
       	 }
       	 else{
       		 if((score > bestScore) && (newMotif.length() == 10)){
       			 bestScore = score;
       			 motif = newMotif;
       			 bestIndex = start + 1;
       		 }
       	 }
       	 start++;
        } // end of while
	    ypResult[0] = motif;
	    ypResult[1] = "" + bestScore;
	    ypResult[2] = "" + bestIndex;
	}  // end of YP Method


	// ***********************************************
	public static void CE(String seq, double[][] scoreMatrix , String[] inrResult, String filePath, int len){


	len = len - 1;
        int start = 0, end = seq.length() - 1, bestIndex = 0;
        double score = 0,bestScore = 0;
        String motif = "",newMotif = "";
        while(start <= end){
       	 score = 0;
       	 newMotif = "";
       	 for(int i = start; i <= (start + len); i++){
       		 if(i < seq.length()){
       			 newMotif+=seq.charAt(i);
       			if(seq.charAt(i) == 'A')
       				score += scoreMatrix[0][i - start];
       			else if(seq.charAt(i) == 'C')
       				score += scoreMatrix[1][i - start];
       			else if(seq.charAt(i) == 'G')
       				score += scoreMatrix[2][i - start];
       			else{
       				if(seq.charAt(i) == 'T')
       			     	score += scoreMatrix[3][i - start];
       			}
       		 } // end of IF for sequence length
       	 } // end of for


       	 if(bestScore == 0){
       		 bestScore = score;
       		 motif = newMotif;
       		 bestIndex = start + 1;
       	 }
       	 else{
       		 if((score > bestScore) && (newMotif.length() == len)){
       			 bestScore = score;
       			 motif = newMotif;
       			 bestIndex = start + 1;
       		 }
       	 }
       	 start++;
        } // end of while
	    inrResult[0] = motif;
	    inrResult[1] = "" + bestScore;
	    inrResult[2] = "" + bestIndex;
	}  // end of INR Method



	// ***********************************************
	public static void NDE(String seq, double[][] scoreMatrix, String[] inrResult, String filePath, int len)
	{

		int start = 0, end = seq.length() - 1, bestIndex = 0;
                len = len - 1;
		double score = 0, bestScore = 0;
		String motif = "", newMotif = "";
		while (start <= end)
		{
			score = 0;
			newMotif = "";
			for (int i = start; i <= (start + len); i++)
			{
				if (i < seq.length())
				{
					newMotif += seq.charAt(i);
					if (seq.charAt(i) == 'A')
						score += scoreMatrix[0][i - start];
					else if (seq.charAt(i) == 'C')
						score += scoreMatrix[1][i - start];
					else if (seq.charAt(i) == 'G')
						score += scoreMatrix[2][i - start];
					else
					{
						if (seq.charAt(i) == 'T')
							score += scoreMatrix[3][i - start];
					}
				} // end of IF for sequence length
			} // end of for


			if (bestScore == 0)
			{
				bestScore = score;
				motif = newMotif;
				bestIndex = start + 1;
			}
			else
			{
				if ((score > bestScore) && (newMotif.length() == len))
				{
					bestScore = score;
					motif = newMotif;
					bestIndex = start + 1;
				}
			}
			start++;
		} // end of while
		inrResult[0] = motif;
		inrResult[1] = "" + bestScore;
		inrResult[2] = "" + bestIndex;
	}  // end of INR Method


	// ***********************************************
	public static void CS(String seq, double[][] scoreMatrix, String[] inrResult, String filePath)
	{

		int start = 147, end = 153, len = 7 - 1, bestIndex = 0;
		double score = 0, bestScore = 0;
		String motif = "", newMotif = "";
		while (start <= end)
		{
			score = 0;
			newMotif = "";
			for (int i = start; i <= (start + len); i++)
			{
				if (i < seq.length())
				{
					newMotif += seq.charAt(i);
					if (seq.charAt(i) == 'A')
						score += scoreMatrix[0][i - start];
					else if (seq.charAt(i) == 'C')
						score += scoreMatrix[1][i - start];
					else if (seq.charAt(i) == 'G')
						score += scoreMatrix[2][i - start];
					else
					{
						if (seq.charAt(i) == 'T')
							score += scoreMatrix[3][i - start];
					}
				} // end of IF for sequence length
			} // end of for


			if (bestScore == 0)
			{
				bestScore = score;
				motif = newMotif;
				bestIndex = start + 1;
			}
			else
			{
				if ((score > bestScore) && (newMotif.length() == len))
				{
					bestScore = score;
					motif = newMotif;
					bestIndex = start + 1;
				}
			}
			start++;
		} // end of while
		inrResult[0] = motif;
		inrResult[1] = "" + bestScore;
		inrResult[2] = "" + bestIndex;
	}  // end of INR Method

	//******************************************************************


    //******************************************************************// Hexamer Scoring
    public static double hexaScore(String useRegion, double[][] hexa, int loc) {
        double score = 0;

        char[] arr = new char[useRegion.length()];
        useRegion.getChars(0, useRegion.length(), arr, 0);
        for (int index = 0; index < (useRegion.length() - 4); index++) {
            String nt = new String();
            nt = nt + getNum(arr[index]) + getNum(arr[index + 1]) + getNum(arr[index + 2]) + getNum(arr[index + 3]) + getNum(arr[index + 4]) ;//+ getNum(arr[index + 5]);
            int i = Integer.parseInt(nt);
            score = score + hexa[i][loc];
        }
        return score;

    } // end of Hexamer scoring

    //******************************************************************
    //******************************************************************
    // get Num
    public static int getNum(char a) {
        if (a == 'A') {
            return 1;
        } else if (a == 'T') {
            return 2;
        } else if (a == 'C') {
            return 3;
        } else {
            if (a == 'G') {
                return 4;
            } else {
                return 0;
            }
        }

    }



    //******************************************************************
    // RoundToFour
    public static double roundTwoFour(String score) {
        double d = Double.parseDouble(score);
        long l = (int) Math.round(d * 10000); // truncates
        d = l / 10000.0;

        return d;

    }


    // *****************************************************
    // Calculate U/A Ratios
    public static double ratioUA(String region){
	int u = 0, a = 0;
        for(int i = 0; i < region.length(); i++){
       		if(region.charAt(i) == 'T')
			u++;
       		if(region.charAt(i) == 'A')
			a++;
        }

        if(a==0)
	  a=1;
	return (u/a);
    }

    // *****************************************************
    // Calculate A/U Ratios
    public static double ratioAU(String region){
	int u = 0, a = 0;
        for(int i = 0; i < region.length(); i++){
       		if(region.charAt(i) == 'T')
			u++;
       		if(region.charAt(i) == 'A')
			a++;
        }

	if(u==0)
	  u = 1;
	return (a/u);
    }


	// ************************************************************ //
	public static double csScore(double[][] csMatrix,String csRegion)
	{
		   double score = 0;
		    int k = 0;

		       for(int j = 0; j < csRegion.length(); j++)
		       {
		         if(csRegion.charAt(k) == 'A')
		           score = score + csMatrix[0][j];
		         else if(csRegion.charAt(k) == 'C')
		           score = score + csMatrix[1][j];
		         else if(csRegion.charAt(k) == 'G')
		           score = score + csMatrix[2][j];
		         else
		           if(csRegion.charAt(k) == 'T')
		             score = score + csMatrix[3][j];
		         k++;
		       }


		     return score;
	}


	  // ************************************************************* //
	  public static double nuecsDis(double[] disMatrix,String[] nueResult)
	  {
	    	  double score = 0;
                  double pos = 41 - (Integer.parseInt(nueResult[2]) + 5);
	    	  if(pos <= 5)
	    		  return disMatrix[0];
	    	  else if(pos <= 10)
	    		  return disMatrix[1];
	    	  else if(pos <= 15)
	    		  return disMatrix[2];
	    	  else if(pos <= 20)
	    		  return disMatrix[3];
	    	  else if(pos <= 25)
	    		  return disMatrix[4];
	    	  else if(pos <= 30)
	    		  return disMatrix[5];
	    	  else if(pos <= 35)
	    		  return disMatrix[6];
			  else
				  return disMatrix[7];
		  			//  else if (pas[0][1] <= 40)
	    		//  return disMatrix[7];
	  }

	  // **********************************************************//
		//NUE distance distribution from CS
		public static double dis(double[] nue,int pos,double score){
		  int diff = 31 - (pos+5);
		  double s = 0;
			//if(score < 1.472)
				//diff = 1;
			//System.out.println("Diff = " + diff);
			if((diff <= 10)&&(diff>=5))
				return nue[0];
			else if((diff <= 15)&&(diff>10))
				return nue[1];
			else if((diff <= 20)&&(diff>15))
				return nue[2];
			else if((diff <= 25)&&(diff>20))
				return nue[3];
			else if((diff <= 30)&&(diff>25))
				return nue[4];
			else{
				if(diff > 30)
					return nue[5];
				else
					return nue[6];
			}
		}

		  // **********************************************************//
		//NUE distance distribution from CS
		public static double disFUE(double[] fue,int pos,double score){

			int diff = 141 - (pos+9);
			//if(Double.parseDouble(fueResult[1]) < 5.156)
				//diff = 1;
			//System.out.println("Diff = " + diff);
			if((diff <= 50)&&(diff>=40))
				return fue[0];
			else if((diff <= 60)&&(diff>50))
				return fue[1];
			else if((diff <= 70)&&(diff>60))
				return fue[2];
			else if((diff <= 80)&&(diff>70))
				return fue[3];
			else if((diff <= 90)&&(diff>80))
				return fue[4];
			else if((diff <= 100)&&(diff>90))
				return fue[5];
			else if((diff <= 110)&&(diff>100))
				return fue[6];
			else if((diff <= 120)&&(diff>110))
				return fue[7];
			else if((diff <= 130)&&(diff>120))
				return fue[8];
			else if((diff <= 140)&&(diff>130))
				return fue[9];
			else
				return fue[10];


		} // end of disFUE



}
