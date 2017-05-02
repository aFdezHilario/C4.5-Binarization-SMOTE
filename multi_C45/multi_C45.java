package multi_C45;

/**
 * <p>Title: Algorithm</p>
 *
 * <p>Description: It contains the implementation of the algorithm</p>
 *
 *
 * <p>Company: KEEL </p>
 *
 * @author Alberto Fernandez
 * @company University of Granada
 * @date 10/11/2012
 * @version 1.0
 */

import java.io.IOException;
import org.core.*;

import C45.C45;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.ArrayList;

public class multi_C45 {

	myDataset train, val, test;
	String outputTr, outputTst, fileRB, majorityClass, binarization, header;
	int nClasses, instancesPerLeaf, n_classifiers, preprocessing;
	float confidence;
	double aprioriClassDistribution[];
	boolean pruned, valid[];
	String fichTrain, method[];
	RuleBase[] ruleBaseTree;
	myDataset[] train_sets;
	ArrayList <C45> classifiers;

	public static int NONE = 0, SMOTE = 1;


	private boolean somethingWrong = false; //to check if everything is correct.

	/**
	 * Default constructor
	 */
	public multi_C45() {
	}

	/**
	 * It reads the data from the input files (training, validation and test) and parse all the parameters
	 * from the parameters array.
	 * @param parameters parseParameters It contains the input files, output files and parameters
	 */
	public multi_C45(parseParameters parameters) {

		train = new myDataset();
		val = new myDataset();
		test = new myDataset();
		fichTrain = parameters.getTrainingInputFile();
		try {
			System.out.println("\nReading the training set: " +
					parameters.getTrainingInputFile());
			train.readClassificationSet(parameters.getTrainingInputFile(), true);
			System.out.println("\nReading the validation set: " +
					parameters.getValidationInputFile());
			val.readClassificationSet(parameters.getValidationInputFile(), false);
			System.out.println("\nReading the test set: " +
					parameters.getTestInputFile());
			test.readClassificationSet(parameters.getTestInputFile(), false);
		}
		catch (IOException e) {
			System.err.println(
					"There was a problem while reading the input data-sets: " +
							e);
			somethingWrong = true;
		}

		//We may check if there are some numerical attributes, because our algorithm may not handle them:
		//somethingWrong = somethingWrong || train.hasRealAttributes();
		//somethingWrong = somethingWrong || train.hasMissingAttributes();

		outputTr = parameters.getTrainingOutputFile();
		outputTst = parameters.getTestOutputFile();

		int param = 0;
		fileRB = parameters.getOutputFile(param++);

		//Now we parse the parameters
		pruned = true;
		confidence = Float.parseFloat(parameters.getParameter(param++));
		instancesPerLeaf = Integer.parseInt(parameters.getParameter(param++));

		binarization = parameters.getParameter(param++);
		classifiers = new ArrayList<C45>();

		String prep = parameters.getParameter(param++);
		if (prep.equalsIgnoreCase("NONE")){
			preprocessing = NONE;
		}else if(prep.equalsIgnoreCase("SMOTE")){
			preprocessing = SMOTE;
		}else{
			System.err.println("Preprocesamiento incorrecto!");
			somethingWrong = true;
		}
		method = new String[2];
		method[0] = "None";
		method[1] = "SMOTE-I";


		//To run in SGE environment
		header = parameters.getTestInputFile();

		String[] aux = null;
		aux = header.split("\\.");
		header = aux[aux.length - 2]; //aux.length-1 is the extension
		aux = header.split("/");
		header = aux[aux.length - 1];    

	}

	/**
	 * It launches the algorithm
	 */
	public void execute() {
		if (somethingWrong) { //We do not execute the program
			System.err.println("An error was found, the data-set has missing values.");
			System.err.println("Aborting the program");
			//We should not use the statement: System.exit(-1);
		}
		else {
			//We do here the algorithm's operations
			long t_ini = System.currentTimeMillis();
			nClasses = train.getnClasses();
			majorityClass = train.claseMasFrecuente();

			n_classifiers = nClasses * (nClasses - 1) / 2;
			if (binarization.equals("OVA")){
				n_classifiers = nClasses;  
			}
			valid = new boolean[n_classifiers];
			ruleBaseTree = new RuleBase[n_classifiers];
			train_sets = new myDataset[n_classifiers];

			aprioriClassDistribution = new double[nClasses];
			for (int i = 0; i < nClasses; i++) {
				aprioriClassDistribution[i] = 1.0 * train.numberInstances(i) /
						train.size();
			}

			if (binarization.equals("OVO")){
				for (int i = 0, x = 0; i < nClasses - 1; i++) {
					for (int j = i + 1; j < nClasses; j++) {
						if (i != j) {
							train_sets[x] = new myDataset(train, i, j);
							//train_sets[x].print();
							x++;
						}
					}
				}
			}else{
				for (int i = 0; i < nClasses; i++) {
					train_sets[i] = new myDataset(train, i);
				}
			}

			int x, y;
			x = 0;
			y = 1;
			for (int i = 0; i < n_classifiers; i++) {
				String cadena = new String("");
				if (binarization.equals("OVO")){
					cadena = train.nombreClase(x)+" vs. "+train.nombreClase(y); 
				}else{
					cadena = train.nombreClase(i)+" vs. REST";
				}
				System.out.println("Classsifier -> "+i+"; "+cadena);
				if (!train_sets[i].empty()) {
					String fileTr = "./datos/"+header+".tra";
					Files.writeFile(fileTr, train_sets[i].printDataSet());
					if (preprocessing != NONE){
						System.out.println("Applying Preprocessing...["+x+"]["+y+"]: "+i);
						try{

							SMOTE smote;
							String fileSMOTE = new String();
							String fileTst = new String();
							fileTst = "./datos/"+header+".tst";

							fileSMOTE += "algorithm = SMOTE\ninputData = \"" + fileTr + "\" \"" + fileTr + "\"\n";
							fileSMOTE += "outputData = \"" + fileTr + "\" \"" + fileTst + "\"\n\n";
							fileSMOTE += "seed = 12345678\n";
							if (train.n_minoritaria() > 5){
								fileSMOTE += "Number of Neighbors = 5\n";
							}else{
								fileSMOTE +=  "Number of Neighbours for SMOTE = "+(train.n_minoritaria()-1)+"\n";	
							}
							fileSMOTE += "Type of SMOTE = BOTH \n";
							fileSMOTE += "Balancing = YES\n";
							fileSMOTE += "Quantity of generated examples = 1\n";
							fileSMOTE += "Distance Function = HVDM\n";
							fileSMOTE += "Type of Interpolation = standard\n";
							fileSMOTE += "Alpha = 0.5\n";
							fileSMOTE += "Mu = 0.5\n";

							//Files.writeFile("./data/"+fileConfName + "_SMOTE_cfg.txt", fileSMOTE);
							//smote = new SMOTE ("./data/"+fileConfName + "_SMOTE_cfg.txt");
							smote = new SMOTE(fileSMOTE);
							smote.run();
							train_sets[i].readClassificationSet(fileTr, false);

							/*
							creaConf(train_sets[i],x);
							Runtime rt = Runtime.getRuntime();
							Process prepro = rt.exec("java -jar ../exe/"+metodo[preprocessing]+".jar ./datos/"+cabecera+".conf");
							//Process smote = rt.exec("java -jar ../exe/SMOTE.jar ./datos/"+cabecera+".conf");
							// any error message?
							StreamGobbler errorGobbler = new
									StreamGobbler(prepro.getErrorStream(), "ERR");

							// any output?
							StreamGobbler outputGobbler = new
									StreamGobbler(prepro.getInputStream(), "OUT");

							// kick them off
							errorGobbler.start();
							outputGobbler.start();
							int valorSalida = prepro.waitFor();
							System.out.println(metodo[preprocessing]+" ejecutado. Pasamos a la fase de aprendizaje " + valorSalida);
							 */
						}
						catch (Throwable t) {
							t.printStackTrace();
						}        	
					}
					valid[i] = true;
					//Files.writeFile(fileTr, train_sets[i].printDataSet());
					//C45 arbol = new C45("./datos/"+cabecera+".tra", pruned, confidence, instancesPerLeaf);
					try{
						C45 tree = new C45(fileTr,pruned, confidence, instancesPerLeaf);
						classifiers.add(tree);
					}catch(Exception e){
						e.printStackTrace();
					}
				}else {
					valid[i] = false;
					System.out.println("No hay ejemplos para la clase "+i);
				}
				y++;
				if (y % nClasses == 0) {
					x++;
					y = x + 1;
				}
			}

			long t_fin = System.currentTimeMillis();
			long t_exec = t_fin - t_ini;
			long hours = t_exec / 3600000;
			long rest = t_exec % 3600000;
			long minutes = rest / 60000;
			rest %= 60000;
			long seconds = rest / 1000;
			rest %= 1000;
			String tiempo = "Execution Time: " + hours + ":" + minutes + ":" + seconds + "." + rest;
			//Finally we should fill the training and test output files
			double accTr = doOutput(this.val, this.outputTr);
			double accTst = doOutput(this.test, this.outputTst);

			writeOutputs(accTr, accTst, tiempo);
		}
	}

	/**
	 * It generates the output file from a given dataset and stores it in a file
	 * @param dataset myDataset input dataset
	 * @param filename String the name of the file
	 * @return the Accuracy of the classifier
	 */
	private double doOutput(myDataset dataset, String filename){
		String output = new String("");
		output = dataset.doHeader(); //we insert the header in the output file
		int [] hits = new int[nClasses];
		//We write the output for each example
		for (int i = 0; i < dataset.getnData(); i++) {
			//for classification:
			String actualClass = dataset.getOutputAsString(i);
			int clase = dataset.getOutputAsInteger(i);
			double [] scores = this.classificationOutput(dataset.getExample(i));
			String prediction = this.getOutputTies(scores);
			output += actualClass + " " + prediction + " "+ scores[this.maxIndex(scores)]+"\n";
			if (actualClass.equalsIgnoreCase(prediction)) {
				hits[clase]++;
			}
		}
		Files.writeFile(filename, output);
		double accAvg = 0;
		int numClases = 0;
		for (int i = 0; i < nClasses; i++){
			try{
				int datos = dataset.numberInstances(i);
				if (datos > 0){
					numClases++;
					double acc = (1.0*hits[i])/datos;
					System.out.print("Cl["+i+"]: "+hits[i]+"/"+datos+"("+acc+")\t");
					accAvg += acc;
				}
			}catch(Exception e){
				System.err.println("No examples for class "+i);
			}
		}
		System.out.println("");
		accAvg /= numClases;
		//accAvg /= nClasses;
		//return (1.0 * hits / dataset.size());
		return (100.0*accAvg);
	}

	/**
	 * It returns the algorithm classification output given an input example
	 * @param example double[] The input example
	 * @return String the output generated by the algorithm
	 */
	private double [] classificationOutput(double[] example) {
		double [] max_scores;
		if (binarization.equals("OVO")){
			double[][] table = computeOVOMatrix(example);
			max_scores = computeClassScoresWV(table);
			//double [] max_scores = computeClassScoresLVPC(tabla);
		}else{
			max_scores = computeClassScoresOVA(example);
		}
		//return this.getOutputMax(max_scores);
		//return this.getOutputTies(max_scores);
		return max_scores;
	}

	private double [][] computeOVOMatrix(double [] example){
		/**
      Here we should include the algorithm directives to generate the
      classification output from the input example
		 */
		double[][] table = new double[nClasses][nClasses];
		int x, y;
		x = 0;
		y = 1;
		//for (C45 tree:classifiers){
		for (int i = 0, j = 0; i < n_classifiers; i++){
			if (valid[i]){
				double votos[] = new double[train.getnClasses()];
				votos = classifiers.get(j++).probs(example);
				table[x][y] = votos[x];
				table[y][x] = votos[y];
			}else{
				table[x][y] = 0;
				table[y][x] = 0;
			}

			y++;
			if (y % nClasses == 0) {
				x++;
				y = x + 1;
			}
		}
		return table;
	}

	private double [] computeClassScoresOVA(double [] example){
		double[] grado_asoc = new double[this.n_classifiers];
		for (int i = 0; i < this.n_classifiers; i++) {
			if (valid[i]) {
				//grado_asoc[i] = baseReglas[i].FRM_OVA(example);
				String clase = new String("?");
				for (int j = 0; (j < ruleBaseTree[i].size()) && (clase.equals("?"));
						j++) {
					if (ruleBaseTree[i].ruleBase.get(j).cubre(example)) {
						clase = ruleBaseTree[i].ruleBase.get(j).clase;
					}
				}
				if (clase.compareTo("positive") == 0){
					grado_asoc[i] = 1;
				}
			}
		}
		return grado_asoc;	  
	}

	private double [] computeClassScoresWV(double [][] tabla){
		double [] max_scores = new double[nClasses];
		for (int i = 0; i < nClasses; i++) {
			double sum_clase = 0.0;
			for (int j = 0; j < nClasses; j++) {
				sum_clase += tabla[i][j];
			}
			max_scores[i] = sum_clase/nClasses;
		}
		return max_scores;
	}

	double[] computeClassScoresLVPC(double[][] tabla) {

		// retrieve the P,C,I matrices
		double[][][] relationSet = relationsForInstance(tabla);
		double[][] preferenceRelation = relationSet[0];
		double[][] conflictRelation = relationSet[1];
		double[][] ignoranceRelation = relationSet[2];

		// calculate the class scores
		double[] classScores = new double[preferenceRelation.length];
		for (int i = 0; i < preferenceRelation.length; i++) {
			for (int j = 0; j < preferenceRelation[i].length; j++) {
				if (i != j) {
					classScores[i] +=
							preferenceRelation[i][j] +
							conflictRelation[i][j] * 0.5 +
							ignoranceRelation[i][j] *
							(aprioriClassDistribution[i] /
									(aprioriClassDistribution[i] +
											aprioriClassDistribution[j]))
							;

				}
			}
			if (Double.isNaN(classScores[i])) {
				classScores[i] = 0;
			}
		}

		if (sum(classScores) == 0) {
			return classScores;
		}

		normalize(classScores);
		return classScores;
	}

	/**
	 * Retrieves the preference, conflict and ignorance matrices for a single instance
	 * @param inst The instance for which the PCI-matrices shall be created
	 * @return [0][][] preference matrix, [1][][] conflict matrix, [2][][] ignorance matrix
	 * @throws Exception
	 */
	public double[][][] relationsForInstance(double[][] tabla) {
		double[][] prefRelation = new double[nClasses][nClasses];
		double[][] conflictRelation = new double[nClasses][nClasses];
		double[][] ignoranceRelation = new double[nClasses][nClasses];

		double s0, s1;

		for (int x = 0; x < nClasses; x++)
			for (int y = 0; y < nClasses; y++)
			{
				s0 = tabla[x][y];
				s1 = tabla[y][x];

				double min = s0 < s1 ? s0 : s1;
				double max = s0 > s1 ? s0 : s1;

				prefRelation[x][y] = s0 - min;
				prefRelation[y][x] = s1 - min;
				conflictRelation[x][y] = min;
				conflictRelation[y][x] = min;
				ignoranceRelation[x][y] = 1 - max;
				ignoranceRelation[y][x] = 1 - max;

				/*   if (s0 == 0 && s1 == 0)
            {
                ignoranceRelation[x][y] = 0;
                ignoranceRelation[y][x] = 0;
            }*/
			}

		return new double[][][] {
			prefRelation, conflictRelation,
			ignoranceRelation};
	}

	public void writeOutputs(double accTr, double accTst, String tiempo) {
		System.out.println("Accuracy in training: " + accTr);
		System.out.println("Accuracy in test: " + accTst);
		System.out.println("Algorithm Finished");
		Fichero.escribeFichero(fileRB,"");
		for (C45 tree:classifiers){
			Files.addToFile(fileRB, tree.toString());
		}
		Files.addToFile(fileRB, "Average Accuracy in training: " + accTr + "\n");
		Files.addToFile(fileRB, "Average Accuracy in test: " + accTst + "\n");
		Files.addToFile(fileRB, tiempo + "\n");
	}

	String getOutputTies(double[] max) {

		/*
		 * Tie-breaking step 1: Find out which classes gain the maximum score
		 */
		double maxValue = max[maxIndex(max)];
		double[] ties = new double[max.length];
		for (int i = 0; i < max.length; i++) {
			if (max[i] == maxValue) {
				ties[i] = aprioriClassDistribution[i];
			}
		}

		max = new double[max.length];
		max[maxIndex(ties)] = 1;

		/*
		 * Tie-breaking step 2: Check whether the tying classes have the same a priori
		 * class probability and count these classes.
		 */
		int tieValues = 0;
		maxValue = ties[maxIndex(ties)];
		for (int i = 0; i < ties.length; i++) {
			if (ties[i] == maxValue) {
				tieValues++;
			}
		}

		/*
		 * Tie-breaking step 3: If the tying classes have the same a priori probabilities,
		 * then use randomization to determine the winner among these classes
		 */
		if (tieValues > 1) {
			tieValues = 0;
			maxValue = ties[maxIndex(ties)];
			int[] stillTying = new int[ties.length];

			for (int i = 0; i < max.length; i++) {
				if (ties[i] == maxValue) {
					stillTying[tieValues] = i;
					tieValues++;
				}
			}
			return train.getOutputValue(stillTying[0]);
		}
		return train.getOutputValue(maxIndex(max));
	}

	String getOutputMax(double[] max) {
		return train.getOutputValue(maxIndex(max));
	}  

	int maxIndex(double[] array) {
		double max = array[0];
		int index = 0;
		for (int i = 1; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
				index = i;
			}
		}
		return index;
	}

	/**
	 * Normalizes the doubles in the array by their sum.
	 *
	 * @param doubles the array of double
	 * @exception IllegalArgumentException if sum is Zero or NaN
	 */
	public static void normalize(double[] doubles) {
		double sum = 0;
		for (int i = 0; i < doubles.length; i++) {
			if (!Double.isNaN(doubles[i]))
				sum += doubles[i];
		}
		normalize(doubles, sum);
	}

	/**
	 * Normalizes the doubles in the array using the given value.
	 *
	 * @param doubles the array of double
	 * @param sum the value by which the doubles are to be normalized
	 * @exception IllegalArgumentException if sum is zero or NaN
	 */
	public static void normalize(double[] doubles, double sum) {

		if (Double.isNaN(sum)) {
			throw new IllegalArgumentException("Can't normalize array. Sum is NaN.");
		}
		if (sum == 0) {
			return;
		}
		for (int i = 0; i < doubles.length; i++) {
			doubles[i] /= sum;
		}
	}

	double sum(double[] array) {
		double sum = 0.0;
		for (int i = 0; i < array.length; i++) {
			sum += array[i];
		}
		return sum;
	}
	  
}
