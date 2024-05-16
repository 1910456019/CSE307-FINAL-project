import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.math.BigDecimal;
import java.math.RoundingMode;

public class PDBParser {
    //
    public static void printMatrix(double[][] matrix) {
        if (matrix == null || matrix.length == 0 || matrix[0].length == 0) {
            System.out.println("The matrix is empty or null");
            return;
        }

        // Find the maximum width of numbers in the matrix for proper alignment
        int maxLength = 0;
        for (double[] row : matrix) {
            for (double num : row) {
                int length = String.format("%.2f", num).length();
                if (length > maxLength) {
                    maxLength = length;
                }
            }
        }

        // Print the matrix with proper alignment
        for (double[] row : matrix) {
            for (double num : row) {
                System.out.printf("%" + (maxLength + 2) + "s", String.format("%.2f", num));
            }
            System.out.println();
        }
    }

    // function that read PDB file and sort the alpha carbons in to hashtable
    public static Map<Integer, double[]> parsePDB(String filePath) throws IOException {
        Map<Integer, double[]> alphaCarbons = new HashMap<>();
        int counter = 1; // Counter for the digit head

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
                    String atomName = line.substring(12, 16).trim();
                    if (atomName.equals("CA")) { // Check if it's an alpha carbon
                        String residueName = line.substring(17, 20).trim();
                        int residueNumber = Integer.parseInt(line.substring(22, 26).trim());
                        double x = Double.parseDouble(line.substring(30, 38).trim());
                        double y = Double.parseDouble(line.substring(38, 46).trim());
                        double z = Double.parseDouble(line.substring(46, 54).trim());
                        Integer residueId = Integer.valueOf(residueNumber); // Use residue number as the digit head

                        if (!alphaCarbons.containsKey(residueId)) {
                            alphaCarbons.put(residueId, new double[] { x, y, z });
                        }
                        counter++; // Increment counter
                    }
                }
            }
        }

        return alphaCarbons;
    }

    public static void visualizeAlphaCarbons(Map<Integer, double[]> alphaCarbons) {
        for (Map.Entry<Integer, double[]> entry : alphaCarbons.entrySet()) {
            Integer residueId = entry.getKey();
            double[] coordinates = entry.getValue();
            System.out.println("Residue ID: " + residueId + ", Coordinates: [" + coordinates[0] + ", " + coordinates[1]
                    + ", " + coordinates[2] + "]");
        }
    }

    // Function to calculate the vector difference between two points in 3D space
    private static double[] calculateDistance(double[] point1, double[] point2) {
        double[] result = new double[3];
        result[0] = point1[0] - point2[0]; // Difference in x-coordinate
        result[1] = point1[1] - point2[1]; // Difference in y-coordinate
        result[2] = point1[2] - point2[2]; // Difference in z-coordinate
        // System.out.println("result: " + result[0] + " " + result[1] + " " +
        // result[2]);
        return result;
    }

    // Function to calculate Euclidean distance squared between two points in 3D
    // space
    private static double calculateDistanceSquared(double[] point1, double[] point2) {
        return Math.pow(point1[0] - point2[0], 2) + Math.pow(point1[1] - point2[1], 2)
                + Math.pow(point1[2] - point2[2], 2);
    }

    private static double calculateLowMatrixValue(double ijrs) {
        int a = 2;
        int b = 100;
        return b / (a + ijrs);
    }

    // compute single low level matrix
    public static double[][] createLowMatrix(Map<Integer, double[]> alphaCarbons1, Map<Integer, double[]> alphaCarbons2,
            int i, int j) {
        // define some value
        // centerPa means picking this specific amino acid in the protein as a center.
        // vectorToCenterPa means the vector from the aminoc acid to the center in same
        // protein
        double[] centerPA = new double[3];
        double[] vectorToCenterPA = new double[3];
        double[] computeVectorPA = new double[3];

        double[] centerPB = new double[3];
        double[] vectorToCenterPB = new double[3];
        double[] computeVectorPB = new double[3];

        // the delta value between two amino acids
        double deltaAminoAcids = 0;

        // Determine the sizes of the sequences
        int proteinALength = alphaCarbons1.size();
        int proteinBLength = alphaCarbons2.size();
        // System.out.println("proteinALength: " + proteinALength);
        // System.out.println("proteinBLength: " + proteinBLength);

        // Initialize the Low level matrix
        double[][] lowMatrix = new double[proteinALength][proteinBLength];

        for (int r = 0; r < proteinALength; r++) {

            for (int s = 0; s < proteinBLength; s++) {
                if (r < i && s >= j || r <= i && s > j) {
                    lowMatrix[r][s] = -1;
                } else if (r > i && s <= j || r >= i && s < j) {
                    lowMatrix[r][s] = -1;
                } else if (r == i && s == j) {
                    // define the vector use to compute the distance between the two amino acids
                    centerPA = alphaCarbons1.get(i + 1);
                    vectorToCenterPA = alphaCarbons1.get(r + 1);

                    centerPB = alphaCarbons2.get(j + 1);
                    vectorToCenterPB = alphaCarbons2.get(s + 1);

                    computeVectorPA = calculateDistance(centerPA, vectorToCenterPA);
                    computeVectorPB = calculateDistance(centerPB, vectorToCenterPB);

                    // compute delta^(i,j)(r,s)
                    deltaAminoAcids = calculateDistanceSquared(computeVectorPA, computeVectorPB);
                    lowMatrix[r][s] = calculateLowMatrixValue(deltaAminoAcids);
                } else {
                    centerPA = alphaCarbons1.get(i + 1);
                    vectorToCenterPA = alphaCarbons1.get(r + 1);

                    centerPB = alphaCarbons2.get(j + 1);
                    vectorToCenterPB = alphaCarbons2.get(s + 1);

                    computeVectorPA = calculateDistance(centerPA, vectorToCenterPA);
                    computeVectorPB = calculateDistance(centerPB, vectorToCenterPB);
                    // System.out.println("computeVectorPA: " + computeVectorPA[0] + " " +
                    // computeVectorPA[1] + " " + computeVectorPA[2]);
                    // System.out.println("computeVectorPB: " + computeVectorPB[0] + " " +
                    // computeVectorPB[1] + " " + computeVectorPB[2]);

                    // compute delta^(i,j)(r,s)
                    deltaAminoAcids = calculateDistanceSquared(computeVectorPA, computeVectorPB);
                    lowMatrix[r][s] = calculateLowMatrixValue(deltaAminoAcids);
                    // System.out.println("deltaAminoAcids: " + deltaAminoAcids);
                    // System.out.println("lowMatrix: " + lowMatrix[r][s]);
                    // System.out.println("lowMatrix: " + lowMatrix[1][1]);
                    // lowMatrix[r][s] = 1;
                }
            }

        }
        // Initialize the mat
        return lowMatrix;
    }

    // compute the High level matrix
    // compute single low level matrix
    public static double[][] createHighMatrix(double[][] lowMatrix, int iIndex, int jIndex) {
        // define high level matrix
        double[][] highMatrix = new double[lowMatrix.length][lowMatrix[0].length];

        int rows = lowMatrix.length;
        int cols = lowMatrix[0].length;
        // define the constant value a,b and gap penalty g.
        int a = 100;
        int b = 2;
        int g = 2;

        // initial base case, which compute the column and row of the center(i,j) index
        // in the high level matrix.
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (j == jIndex) {
                    highMatrix[i][j] = (double) a / b - g * Math.abs(i - iIndex);
                } else if (i == iIndex) {
                    highMatrix[i][j] = (double) a / b - g * Math.abs(j - jIndex);
                } else {
                    highMatrix[i][j] = -1;
                }
            }
        }

        // Recurrence relation for the bottom-right submatrix (r > i, s > j) &
        // bottom-left submatrix (r < i, s < j)
        for (int r = 0; r < rows; r++) {
            for (int s = 0; s < cols; s++) {
                // Recurrence relation for the bottom-right submatrix (r > i, s > j)
                if (r > iIndex && s > jIndex) {
                    highMatrix[r][s] = Math.max(
                            highMatrix[r - 1][s - 1] + lowMatrix[r][s],
                            Math.max(
                                    highMatrix[r - 1][s] - g,
                                    highMatrix[r][s - 1] - g));
                }

                // // Recurrence relation for the bottom-left submatrix (r < i, s < j)
                // if (r < iIndex && s < jIndex) {
                // highMatrix[r][s] = Math.max(
                // highMatrix[r + 1][s + 1] + lowMatrix[r][s],
                // Math.max(
                // highMatrix[r + 1][s] - g,
                // highMatrix[r][s + 1] - g
                // )
                // );

                // }
            }
        }

        // Recurrence relation for the top-left submatrix (r < i, s < j)
        for (int r = rows; r >= 0; r--) {
            for (int s = cols; s >= 0; s--) {
                // Recurrence relation for the bottom-left submatrix (r < i, s < j)
                if (r < iIndex && s < jIndex) {
                    highMatrix[r][s] = Math.max(
                            highMatrix[r + 1][s + 1] + lowMatrix[r][s],
                            Math.max(
                                    highMatrix[r + 1][s] - g,
                                    highMatrix[r][s + 1] - g));

                }

            }
        }
        return highMatrix;
    }

    // implemantation of the Smith-waterman algorithm
    public static double[][] updateMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        // define matrix
        double[][] updateMatrix = new double[rows][cols];
        // initialize the matrix
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {

                updateMatrix[i][j] = 0;

            }
        }
        return updateMatrix;
    }

    public static double[] findAlignment(double[][] matrix, int iIndex, int jIndex) {
        // define values
        List<Double> tempAlignmentPath = new ArrayList<>();
        List<Double> temp2AlignmentPath = new ArrayList<>();
        int rows = matrix.length;
        int cols = matrix[0].length;

        int i = iIndex;
        int j = jIndex;

        while (i >= 0 && j >= 0 && matrix[i][j] > 0) {
            tempAlignmentPath.add(matrix[i][j]);
            if (i > 0 && j > 0 && matrix[i - 1][j - 1] >= matrix[i - 1][j]
                    && matrix[i - 1][j - 1] >= matrix[i][j - 1]) {
                i--;
                j--;
            } else if (i > 0 && (j == 0 || matrix[i - 1][j] >= matrix[i][j - 1])) {
                i--;
            } else if (j > 0) {
                j--;
            } else {
                break;
            }
        }

        i = iIndex;
        j = jIndex;
        // Alignment starting from (iIndex, jIndex) and moving to the bottom-right
        while (i < rows && j < cols && matrix[i][j] > 0) {
            if (!temp2AlignmentPath.contains(matrix[i][j])) {
                temp2AlignmentPath.add(matrix[i][j]);
            }
            if (i < rows - 1 && j < cols - 1 && matrix[i + 1][j + 1] >= matrix[i + 1][j]
                    && matrix[i + 1][j + 1] >= matrix[i][j + 1]) {
                i++;
                j++;
            } else if (i < rows - 1 && (j == cols - 1 || matrix[i + 1][j] >= matrix[i][j + 1])) {
                i++;
            } else if (j < cols - 1) {
                j++;
            } else {
                break;
            }
        }
        temp2AlignmentPath.remove(0);
        Collections.reverse(temp2AlignmentPath);

        // Combine the two alignment paths
        temp2AlignmentPath.addAll(tempAlignmentPath);

        // Convert the temporary list to an array
        double[] alignmentPath = new double[temp2AlignmentPath.size()];
        for (int k = 0; k < temp2AlignmentPath.size(); k++) {
            alignmentPath[k] = temp2AlignmentPath.get(temp2AlignmentPath.size() - k - 1);
        }
        return alignmentPath;
    }

    // find the similarity score for the matrix(two max value in the matrix).
    public static double findSimilarityScore(double[] alignedSequence) {
        // Initialize the two largest numbers to negative infinity
        double max1 = Double.NEGATIVE_INFINITY;
        double max2 = Double.NEGATIVE_INFINITY;

        for (double num : alignedSequence) {
            if (num > max1) {
                // Update max2 before changing max1
                max2 = max1;
                max1 = num;
            } else if (num > max2) {
                max2 = num;
            }
        }

        return max1 + max2;
    }

    // compute the threshold value
    // Method to calculate the mean
    public static double calculateMean(double[] scores) {
        double sum = 0;
        for (double score : scores) {
            sum += score;
        }
        return sum / scores.length;
    }

    // Method to calculate the standard deviation
    public static double calculateStandardDeviation(double[] scores, double mean) {
        double sum = 0;
        for (double score : scores) {
            sum += Math.pow(score - mean, 2);
        }
        return Math.sqrt(sum / scores.length);
    }

    public static double[][] updateScoreMatrix(double[][] matrix, double[][] Smatrix, int iIndex, int jIndex) {
        int rows = matrix.length;
        int cols = matrix[0].length;

        // Ensure the final score matrix is initialized to zeros
        // for (int i = 0; i < rows; i++) {
        // for (int j = 0; j < cols; j++) {
        // finalScoreMatrix[i][j] = 0.0;
        // }
        // }

        int i = iIndex;
        int j = jIndex;

        while (i >= 0 && j >= 0 && matrix[i][j] > 0) {
            Smatrix[i][j] = matrix[i][j];
            if (i > 0 && j > 0 && matrix[i - 1][j - 1] >= matrix[i - 1][j]
                    && matrix[i - 1][j - 1] >= matrix[i][j - 1]) {
                i--;
                j--;
            } else if (i > 0 && (j == 0 || matrix[i - 1][j] >= matrix[i][j - 1])) {
                i--;
            } else if (j > 0) {
                j--;
            } else {
                break;
            }
        }

        i = iIndex;
        j = jIndex;
        // Alignment starting from (iIndex, jIndex) and moving to the bottom-right
        while (i < rows && j < cols && matrix[i][j] > 0) {
            Smatrix[i][j] = matrix[i][j];
            if (i < rows - 1 && j < cols - 1 && matrix[i + 1][j + 1] >= matrix[i + 1][j]
                    && matrix[i + 1][j + 1] >= matrix[i][j + 1]) {
                i++;
                j++;
            } else if (i < rows - 1 && (j == cols - 1 || matrix[i + 1][j] >= matrix[i][j + 1])) {
                i++;
            } else if (j < cols - 1) {
                j++;
            } else {
                break;
            }
        }
        return Smatrix;
    }

    // Function to align two sequences using the Smith-Waterman algorithm
    public static void alignSequences(double[][] scoreMatrix, String seq1, String seq2) {
        int rows = seq1.length() + 1;
        int cols = seq2.length() + 1;

        double[][] H = new double[rows][cols];
        int[][] traceback = new int[rows][cols];

        double maxScore = 0;
        int maxI = 0;
        int maxJ = 0;

        for (int i = 1; i < rows; i++) {
            for (int j = 1; j < cols; j++) {
                double match = H[i - 1][j - 1] + scoreMatrix[i - 1][j - 1];
                double delete = H[i - 1][j] + 10;
                double insert = H[i][j - 1] + 10;
                H[i][j] = Math.max(0, Math.max(match, Math.max(delete, insert)));

                if (H[i][j] == match) {
                    traceback[i][j] = 1; // Diagonal
                } else if (H[i][j] == delete) {
                    traceback[i][j] = 2; // Up
                } else if (H[i][j] == insert) {
                    traceback[i][j] = 3; // Left
                }

                if (H[i][j] > maxScore) {
                    maxScore = H[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        StringBuilder align1 = new StringBuilder();
        StringBuilder align2 = new StringBuilder();

        int i = maxI;
        int j = maxJ;

        while (i > 0 && j > 0 && H[i][j] > 0) {
            if (traceback[i][j] == 1) {
                align1.append(seq1.charAt(i - 1));
                align2.append(seq2.charAt(j - 1));
                i--;
                j--;
            } else if (traceback[i][j] == 2) {
                align1.append(seq1.charAt(i - 1));
                align2.append('-');
                i--;
            } else if (traceback[i][j] == 3) {
                align1.append('-');
                align2.append(seq2.charAt(j - 1));
                j--;
            }
        }

        System.out.println("Aligned Sequences:");
        System.out.println(align1.reverse().toString());
        System.out.println(align2.reverse().toString());
    }

    // public static void aliMatrix(double[][] matrix, int iIndex, int jIndex) {
    // int rows = matrix.length;
    // int cols = matrix[0].length;

    // int i = iIndex;
    // int j = jIndex;

    // while (i >= 0 && j >= 0 && matrix[i][j] > 0) {
    // Smatrix[i][j] = matrix[i][j];
    // if (i > 0 && j > 0 && matrix[i-1][j-1] >= matrix[i-1][j] && matrix[i-1][j-1]
    // >= matrix[i][j-1]) {
    // i--;
    // j--;
    // } else if (i > 0 && (j == 0 || matrix[i-1][j] >= matrix[i][j-1])) {
    // i--;
    // } else if (j > 0) {
    // j--;
    // } else {
    // break;
    // }
    // }

    // i = iIndex;
    // j = jIndex;
    // // Alignment starting from (iIndex, jIndex) and moving to the bottom-right
    // while (i < rows && j < cols && matrix[i][j] > 0) {
    // Smatrix[i][j] = matrix[i][j];
    // if (i < rows - 1 && j < cols - 1 && matrix[i + 1][j + 1] >= matrix[i + 1][j]
    // && matrix[i + 1][j + 1] >= matrix[i][j + 1]) {
    // i++;
    // j++;
    // } else if (i < rows - 1 && (j == cols - 1 || matrix[i + 1][j] >= matrix[i][j
    // + 1])) {
    // i++;
    // } else if (j < cols - 1) {
    // j++;
    // } else {
    // break;
    // }
    // }
    // }
    public static int[] findMaxValueIndex(double[][] matrix) {
        int[] maxIndex = new int[2];
        double maxValue = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > maxValue) {
                    maxValue = matrix[i][j];
                    maxIndex[0] = i;
                    maxIndex[1] = j;
                }
            }
        }
        return maxIndex;
    }

    public static String[] needlemanWunschTraceback(double[][] scoreMatrix, String seq1, String seq2) {
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();

        int i = seq2.length();
        int j = seq1.length();

        // // Special case to handle the last digit
        // alignedSeq1.insert(0, seq1.charAt(j));
        // alignedSeq2.insert(0, seq2.charAt(i));
        // i--;
        // j--;

        while (i > 0 || j > 0) {
            if (i > 1 && j > 1 && (scoreMatrix[i - 1][j - 1] == scoreMatrix[i - 2][j - 2])) {
                alignedSeq1.insert(0, seq1.charAt(j - 1));
                alignedSeq2.insert(0, seq2.charAt(i - 1));
                i--;
                j--;
            } else if (i > 1 && j > 1 && i < seq2.length() && j < seq1.length()
                    && (scoreMatrix[i - 1][j - 1] == scoreMatrix[i][j])) {
                alignedSeq1.insert(0, seq1.charAt(j - 1));
                alignedSeq2.insert(0, seq2.charAt(i - 1));
                i--;
                j--;
            } else if (i > 0 && (j == 0 || scoreMatrix[i - 1][j - 1] == scoreMatrix[i - 2][j - 1])) {
                alignedSeq1.insert(0, '-');
                alignedSeq2.insert(0, seq2.charAt(i - 1));
                i--;
            } else {
                alignedSeq1.insert(0, seq1.charAt(j - 1));
                alignedSeq2.insert(0, '-');
                j--;
            }
        }

        // while (i > 0 || j > 0) {
        // if (i > 1 && j > 1 && (scoreMatrix[i-1][j-1] != 0 && scoreMatrix[i - 2][j -
        // 2] != 0)) {
        // alignedSeq1.insert(0, seq1.charAt(j - 1));
        // alignedSeq2.insert(0, seq2.charAt(i - 1));
        // i--;
        // j--;
        // } else if (i > 1 && j > 1 && i < seq2.length() - 1 && j < seq1.length() - 1
        // && (scoreMatrix[i-1][j-1] != 0 && scoreMatrix[i][j] != 0)) {
        // alignedSeq1.insert(0, seq1.charAt(j - 1));
        // alignedSeq2.insert(0, seq2.charAt(i - 1));
        // i--;
        // j--;
        // } else if (i > 1 && j > 1 && scoreMatrix[i-1][j-1] != 0 && scoreMatrix[i -
        // 1][j-2] != 0) {
        // alignedSeq1.insert(0, '-');
        // alignedSeq2.insert(0, seq2.charAt(i - 1));
        // i--;
        // } else if (i > 1 && j > 1) {
        // alignedSeq1.insert(0, seq1.charAt(j - 1));
        // alignedSeq2.insert(0, '-');
        // j--;
        // }
        // }

        return new String[] { alignedSeq1.toString(), alignedSeq2.toString() };
    }

    public static void printAlignedSequences(String[] alignedSequences) {
        if (alignedSequences == null || alignedSequences.length != 2) {
            System.out.println("Invalid aligned sequences array.");
            return;
        }
        System.out.println("Aligned Sequence 1: " + alignedSequences[0]);
        System.out.println("Aligned Sequence 2: " + alignedSequences[1]);
    }

    public static double sumMaxScores(double[][] matrix) {
        double max1 = Double.NEGATIVE_INFINITY;
        double max2 = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > max1) {
                    max2 = max1;
                    max1 = matrix[i][j];
                } else if (matrix[i][j] > max2) {
                    max2 = matrix[i][j];
                }
            }
        }

        return max1 + max2;
    }

    public static void main(String[] args) {
        String pdbFilePath1 = "1l2y.pdb";
        String pdbFilePath2 = "2m7d.pdb";

        String seq2 = "NLYIQWLKDGGPSSGRPPPS";
        String seq1 = "DAYAQWLADAGWASARPPPS";

        try {
            // define the alpha carbons for each protein
            Map<Integer, double[]> alphaCarbons1 = parsePDB(pdbFilePath1);
            Map<Integer, double[]> alphaCarbons2 = parsePDB(pdbFilePath2);

            // define matrix length and similarity score
            int proteinALength = alphaCarbons1.size();
            int proteinBLength = alphaCarbons2.size();

            double[] similarityScores = new double[proteinALength * proteinBLength];

            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            for (int i = 0; i < proteinALength; i++) {
                for (int j = 0; j < proteinBLength; j++) {
                    // for each center pair (i,j) compute the low level matrix first
                    double[][] lowMatrix = createLowMatrix(alphaCarbons1, alphaCarbons2, i, j);
                    // then compute the high level matrix by using low level matrix
                    double[][] highMatrix = createHighMatrix(lowMatrix, i, j);
                    // then Aligment the high level matrix to find the socre you want to use to
                    // update the matrix
                    double[] alignmentResult = findAlignment(highMatrix, i, j);
                    // before you know you want to update the matrix or not, you want to know the
                    // similarity score
                    double similarityScore = findSimilarityScore(alignmentResult);
                    // if the similarity score is greater above than the threshold, update the
                    // matrix
                    similarityScores[i * proteinBLength + j] = similarityScore;

                }
            }

            // Calculate the mean
            double mean = calculateMean(similarityScores);
            // System.out.println("Mean: " + mean);

            // Calculate the standard deviation
            double standardDeviation = calculateStandardDeviation(similarityScores, mean);
            // System.out.println("Standard Deviation: " + standardDeviation);

            // Compute the threshold
            double threshold = mean + standardDeviation;
            // System.out.println("Threshold (Mean + 1*Standard Deviation): " + threshold);

            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            // ----------------------------------compute the threshold for simialirity
            // score----------------------------------
            // ----------------------------------END----------------------------------

            // ----------------------------------main
            // logic----------------------------------
            // define the final matrix
            double[][] finalMatrix = new double[proteinALength][proteinBLength];
            // Ensure the final score matrix is initialized to zeros
            for (int i = 0; i < proteinALength; i++) {
                for (int j = 0; j < proteinBLength; j++) {
                    finalMatrix[i][j] = 0.0;
                }
            }

            // Run the double dynamic to get the final score matrix
            for (int i = 0; i < proteinALength; i++) {
                for (int j = 0; j < proteinBLength; j++) {
                    // for each center pair (i,j) compute the low level matrix first
                    double[][] lowMatrix = createLowMatrix(alphaCarbons1, alphaCarbons2, i, j);
                    // then compute the high level matrix by using low level matrix
                    double[][] highMatrix = createHighMatrix(lowMatrix, i, j);
                    // then Aligment the high level matrix to find the socre you want to use to
                    // update the matrix
                    double[] alignmentResult = findAlignment(highMatrix, i, j);
                    // before you know you want to update the matrix or not, you want to know the
                    // similarity score
                    double similarityScore = findSimilarityScore(alignmentResult);
                    // if the similarity score is greater above than the threshold, update the
                    // matrix
                    if (similarityScore >= threshold) {
                        finalMatrix = updateScoreMatrix(highMatrix, finalMatrix, i, j);
                    }
                }
            }
            // printMatrix(finalMatrix);
            System.out.println();
            int[] index = findMaxValueIndex(finalMatrix);
            double[][] NonMatrix = new double[proteinALength][proteinBLength];

            double[][] FM = updateScoreMatrix(finalMatrix, NonMatrix, index[0], index[1]);
            // printMatrix(FM);
            // clean matrix
            for (int i = 0; i < proteinALength - 1; i++) {
                for (int j = 0; j < proteinBLength - 1; j++) {
                    if (FM[i][j] != 0.0 && FM[i + 1][j + 1] != 0.0) {
                        FM[i][j + 1] = 0.0;
                        FM[i + 1][j] = 0.0;
                    }
                }
            }
            // System.out.println();
            printMatrix(FM);

            String[] output = needlemanWunschTraceback(FM, seq1, seq2);
            printAlignedSequences(output);
            double similarityScore = sumMaxScores(FM);
            System.out.println("Similarity Score for aligment: " + similarityScore);
            // ----------------------------------main
            // logic(END)----------------------------------

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
