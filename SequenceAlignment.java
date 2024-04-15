import java.io.*;
import java.util.List;
import java.util.ArrayList;


public class SequenceAlignment {
    public static String whatMatrix = "";
    public static String backSeq1 = "";
    public static String backSeq2 = "";
    public static int icounter;
    public static int jcounter;


    public static void main(String[] args) throws IOException {
        // Read in the reads from a Fasta file
        List<String> sequence = readReadsFromFastaFile("shortSeq.fasta");


        // split the sequence to two string.
        String seq1 = sequence.get(0);
        String seq2 = sequence.get(1);


        //set the length for seq1 and seq2 as a i,j global counter
        icounter = seq2.length();
        jcounter = seq1.length();


        // define Name for print matrix
        String nameM = "M -- Matrix";
        String nameX = "X -- Matrix";
        String nameY = "Y -- Matrix";


        // define matrix
        double[][] matrixM = new double[seq2.length() + 1][seq1.length() + 1];
        double[][] matrixX = new double[seq2.length() + 1][seq1.length() + 1];
        double[][] matrixY = new double[seq2.length() + 1][seq1.length() + 1];


        // build matrix
        buildMatrix(matrixM, matrixX, matrixY, seq1, seq2);


        // print matrix
        printMatrix(matrixM, nameM, seq1, seq2);
        printMatrix(matrixX, nameX, seq1, seq2);
        printMatrix(matrixY, nameY, seq1, seq2);


        backtrack(matrixM,matrixX,matrixY,seq1,seq2);


        System.out.println(backSeq1);
        System.out.println(backSeq2);
        // Write the sequence to a Fasta file
        // writeContigsToFastaFile(sequence, "output.fasta");
    }


    public static void buildMatrix(double[][] m, double[][] x, double[][] y, String seq1, String seq2) {
        buildBaseMatrix(m, x, y);


        // socring all the matrix
        for (int i = 1; i < m.length; i++) {
            for (int j = 1; j < m[0].length; j++) {
                // score the X matrix
                initX(m, x, y, i, j);


                // score the X matrix
                initY(m, x, y, i, j);


                // score the X matrix
                initM(m, x, y, i, j, seq1, seq2);
            }
        }
    }


    // inpelemnt the base for M,X,Y
    public static void buildBaseMatrix(double[][] m, double[][] x, double[][] y) {
        // initial the base in M matrix
        initBaseM(m);
        // initial the base in X matrix
        initBaseX(x);
        // initial the base in Y matrix
        initBaseY(y);
    }


    // backtrack distinguish
    public static void backtrack(double[][] m, double[][] x, double[][] y, String seq1, String seq2) {
        while (icounter > 0 && jcounter > 0) {
            System.out.println("in the while ---" + whatMatrix);
            if (whatMatrix.equals("")) {
                backtrackAtM(icounter,jcounter,m,x,y,seq1,seq2);
            }
            if (whatMatrix.equals("M")) {
                backtrackAtM(icounter,jcounter,m,x,y,seq1,seq2);
                // System.out.println(icounter);
                // System.out.println(jcounter);
                // System.out.println(backSeq1);
                // System.out.println(backSeq2);
            }
            if (whatMatrix.equals("X")) {
                backtrackAtM(icounter,jcounter,m,x,y,seq1,seq2);
                // System.out.println(icounter);
                // System.out.println(jcounter);
                // System.out.println(backSeq1);
                // System.out.println(backSeq2);
            }
            if (whatMatrix.equals("Y")) {
                backtrackAtM(icounter,jcounter,m,x,y,seq1,seq2);
                // System.out.println(icounter);
                // System.out.println(jcounter);
                // System.out.println(backSeq1);
                // System.out.println(backSeq2);
            }
        }
        String tempStr1 = reversefunction(backSeq1);
        String tempStr2 = reversefunction(backSeq2);
        backSeq1 = tempStr1;
        backSeq2 = tempStr2;
    }


    // backtrack at M matrix
    public static void backtrackAtM(int i, int j, double[][] m, double[][] x, double[][] y, String seq1, String seq2) {
        double score1, score2, score3, matchVal, maxVal;
        String tempStr1;
        String tempStr2;
        System.out.println("I'm in M--matrix");
        System.out.println("i index --" + i);
        System.out.println("j index --" + j);


        if (i == 0 && j == 0) {
            System.out.print("backtrack over");
            return;
        }


        if (i > 0 && j > 0) {
            // find three score
            score1 = m[i-1][j-1];
            score2 = x[i][j];
            score3 = y[i][j];
            if (i == 1 && j == 1) {
                matchVal = 1;
            } else {
                matchVal = matches(i - 1, j - 1, seq1, seq2);
            }


            System.out.println(score1+matchVal);
            System.out.println(score2);
            System.out.println(score3);


            maxVal = Math.max(Math.max(score1 + matchVal, score2), score3);
            System.out.print("max score --" + maxVal);
            System.out.println();


            if (maxVal == x[i][j]) {
                System.out.println("I'm here");


                whatMatrix = "X";


                System.out.println(whatMatrix);
                System.out.println();


                tempStr1 = backSeq1 + seq1.charAt(j - 1);
                tempStr2 = backSeq2 + "-";
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                // System.out.println(backSeq1);
                // System.out.println(backSeq2);
                jcounter -= 1;
            } else if (maxVal == y[i][j]) {
                whatMatrix = "Y";
                tempStr1 = backSeq1 + "-";
                tempStr2 = backSeq2 + seq2.charAt(i - 1);
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                // System.out.println(backSeq1);
                // System.out.println(backSeq2);
                icounter -= 1;
            } else {
                whatMatrix = "M";
                tempStr1 = backSeq1 + seq1.charAt(j - 1);
                tempStr2 = backSeq2 + seq2.charAt(i - 1);
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                // System.out.println(backSeq1);
                // System.out.println(backSeq2);
                icounter -= 1;
                jcounter -= 1;
            }


        }
    }


    // backtrack at X matrix
    public static void backtrackAtX(int i, int j, double[][] m, double[][] x, double[][] y, String seq1, String seq2) {
        double score1, score2, score3, maxVal;
        String tempStr1;
        String tempStr2;
        double S = -10;
        double E = -0.5;


        System.out.println("I'm in X--matrix");
        System.out.println("i index --" + i);
        System.out.println("j index --" + j);


        if (i == 0 && j == 0) {
            System.out.print("backtrack over");
            return;
        }
        if (i > 0 && j > 0) {
            // find three score
            score1 = S + E + m[i][j];
            score2 = E + x[i][j];
            score3 = S + E + y[i][j];


            // pick the max score
            maxVal = Math.max(Math.max(score1, score2), score3);
            System.out.print("max score --" + maxVal);


            if (maxVal == score2) {
                whatMatrix = "X";
                tempStr1 = backSeq1 + seq1.charAt(j - 1);
                tempStr2 = backSeq2 + "-";
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                System.out.println(backSeq1);
                System.out.println(backSeq2);
                jcounter -= 1;
            } else if (maxVal == score3) {
                whatMatrix = "Y";
                tempStr1 = backSeq1 + "-";
                tempStr2 = backSeq2 + seq2.charAt(i - 1);
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                System.out.println(backSeq1);
                System.out.println(backSeq2);
                jcounter -= 1;
            } else {
                whatMatrix = "M";
                tempStr1 = backSeq1 + seq1.charAt(j - 1);
                tempStr2 = backSeq2 + seq2.charAt(i - 1);
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                System.out.println(backSeq1);
                System.out.println(backSeq2);
                jcounter -= 1;
            }
        }
    }


    // backtrack at Y matrix
    public static void backtrackAtY(int i, int j, double[][] m, double[][] x, double[][] y, String seq1, String seq2) {
        double score1, score2, score3, maxVal;
        String tempStr1;
        String tempStr2;
        double S = -10;
        double E = -0.5;
       
        System.out.println("I'm in Y--matrix");
        System.out.println("i index --" + i);
        System.out.println("j index --" + j);
       
        if (i == 0 && j == 0) {
            System.out.print("backtrack over");
            return;
        }
        if (i > 0 && j > 0) {
            // find three score
            score1 = S + E + m[i][j];
            score2 = S + E + x[i][j];
            score3 = E + y[i][j];


            // pick the max score
            maxVal = Math.max(Math.max(score1, score2), score3);
           
            if (maxVal == score2) {
                whatMatrix = "X";
                tempStr1 = backSeq1 + seq1.charAt(j - 1);
                tempStr2 = backSeq2 + "-";
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                icounter -= 1;
            } else if (maxVal == score3) {
                whatMatrix = "Y";
                tempStr1 = backSeq1 + "-";
                tempStr2 = backSeq2 + seq2.charAt(i - 1);
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                icounter -= 1;
            } else {
                whatMatrix = "M";
                tempStr1 = backSeq1 + seq1.charAt(j - 1);
                tempStr2 = backSeq2 + seq2.charAt(i - 1);
                backSeq1 = tempStr1;
                backSeq2 = tempStr2;
                icounter -= 1;
            }
        }
    }


    // mached
    public static int matches(int i, int j, String seq1, String seq2) {
        // System.out.print(seq2.charAt(i-1));
        // System.out.print(seq1.charAt(j-1));
        // System.out.println();
        if (seq2.charAt(i - 1) == seq1.charAt(j - 1)) {
            return 1;
        } else {
            return -4;
        }
    }


    // implement base row and col for matrix M
    public static void initBaseM(double[][] m) {
        int numRows = m.length;
        int numCols = m[0].length;
        // fill the first element as 0
        m[0][0] = 0;
        // fill in rest of first row
        for (int j = 1; j < numCols; j++) {
            m[0][j] = Double.NEGATIVE_INFINITY;
        }
        for (int i = 1; i < numRows; i++) {
            m[i][0] = Double.NEGATIVE_INFINITY;
        }
    }


    // implement base row and col for matrix Y
    public static void initBaseY(double[][] y) {
        int numRows = y.length;
        int numCols = y[0].length;
        double S = -10;
        double E = -0.5;
        // fill the first element as 0
        y[0][0] = 0;
        // fill in rest of first row
        for (int j = 1; j < numCols; j++) {
            y[0][j] = Double.NEGATIVE_INFINITY;
        }
        for (int i = 1; i < numRows; i++) {
            y[i][0] = S + E;
            E += -0.5;
        }
    }


    // implement base row and col for matrix X
    public static void initBaseX(double[][] x) {
        int numRows = x.length;
        int numCols = x[0].length;
        double S = -10;
        double E = -0.5;
        // fill the first element as 0
        x[0][0] = 0;
        // fill in rest of first row
        for (int j = 1; j < numCols; j++) {
            x[0][j] = S + E;
            E += -0.5;
        }
        for (int i = 1; i < numRows; i++) {
            x[i][0] = Double.NEGATIVE_INFINITY;
        }
    }


    // implement matrix M
    public static void initM(double[][] m, double[][] x, double[][] y, int i, int j, String seq1, String seq2) {
        // define varable need to used
        double score1, score2, score3;
        double maxVal = 0;


        // find three score
        score1 = m[i - 1][j - 1];
        score2 = x[i][j];
        score3 = y[i][j];


        // pick the max score
        maxVal = Math.max(Math.max(score1, score2), score3);
        m[i][j] = maxVal + matches(i, j, seq1, seq2);
    }


    // implement matrix X
    public static void initX(double[][] m, double[][] x, double[][] y, int i, int j) {
        double score1, score2, score3;
        double S = -10;
        double E = -0.5;
        // m[0][0] = Double.NEGATIVE_INFINITY;
        // x[0][0] = Double.NEGATIVE_INFINITY;
        // y[0][0] = Double.NEGATIVE_INFINITY;
        score1 = S + E + m[i][j - 1];
        score2 = E + x[i][j - 1];
        score3 = S + E + y[i][j - 1];
        // System.out.println("cordinate value --" + m[i][j-1] + "--j cor --" + j);
        // System.out.println("cordinate value --" + x[i][j-1] + "--j cor --" + j );
        // System.out.println("cordinate value --" + y[i][j-1] + "--j cor --" + j);


        x[i][j] = Math.max(Math.max(score1, score2), score3);
        // m[0][0] = 0;
        // x[0][0] = 0;
        // y[0][0] = 0;


        // if (Double.isInfinite(m[0][j-1])) {
        // score1 = Double.NEGATIVE_INFINITY;
        // } else {
        // score1 = S + E + m[0][j-1];
        // }
        // if (Double.isInfinite(x[0][j-1])) {
        // score2 = Double.NEGATIVE_INFINITY;
        // } else {
        // score2 = E + x[0][j-1];
        // }
        // if (Double.isInfinite(y[0][j-1])) {
        // score3 = Double.NEGATIVE_INFINITY;
        // } else {
        // score3 = S + E + y[0][j-1];
        // }
        // System.out.println("score-1" + score1);
        // System.out.println("score-2" + score2);
        // System.out.println("score-3" + score3);
        // x[i][j] = Math.max(Math.max(score1, score2), score3);
    }


    // implement matrix Y
    public static void initY(double[][] m, double[][] x, double[][] y, int i, int j) {
        double score1, score2, score3;
        double S = -10;
        double E = -0.5;
        // find each score
        score1 = S + E + m[i - 1][j];
        score2 = S + E + x[i - 1][j];
        score3 = E + y[i - 1][j];
        // choose max score
        y[i][j] = Math.max(Math.max(score1, score2), score3);
    }


    private static List<String> readReadsFromFastaFile(String filename) throws IOException {
        List<String> reads = new ArrayList<>();


        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            StringBuilder sequence = new StringBuilder();
            String line;


            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    // reached a new reads, so save the previous one
                    if (sequence.length() > 0) {
                        reads.add(sequence.toString());
                        sequence.setLength(0);
                    }
                } else {
                    // Add the current line to the current reads
                    sequence.append(line.trim());
                }
            }
            // Save the last reads
            if (sequence.length() > 0) {
                reads.add(sequence.toString());
            }
        }
        return reads;
    }


    private static void writeContigsToFastaFile(List<String> contigs, String filename) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            int count = 0;
            for (int i = 0; i < contigs.size(); i++) {
                String header = ">seq_" + (i + 1);
                String sequence = contigs.get(i);
                // Write the header and sequence to the Fasta file
                writer.write(header);
                writer.newLine();
                writer.write(sequence);
                writer.newLine();
                count++;
            }
            // System.out.println(count);
        }
    }


    public static void printMatrix(double[][] matrix, String name, String seq1, String seq2) {
        System.out.println(name + ":");
        // Print the sequence letters for seq1 at the left of the matrix
        System.out.print("\t\t");
        for (int i = 0; i < seq1.length(); i++) {
            System.out.printf("%-7s", seq1.charAt(i) + " ");
        }
        System.out.println();


        for (int i = 0; i < matrix.length; i++) {
            if (i > 0) {
                System.out.print(seq2.charAt(i - 1) + "\t");
            } else {
                System.out.print("\t");
            }
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.printf("%-7.1f", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }


    public static String reversefunction(String str) {
        StringBuilder sb = new StringBuilder(str);
        return sb.reverse().toString();
    }
}
