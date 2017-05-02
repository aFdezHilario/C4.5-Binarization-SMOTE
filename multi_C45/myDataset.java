package multi_C45;

/**
 * <p>Title: Dataset</p>
 *
 * <p>Description: It contains the methods to read a Classification/Regression Dataset</p>
 *
 *
 * <p>Company: KEEL </p>
 *
 * @author Alberto Fern�ndez
 * @version 1.0
 */

import java.io.IOException;

import keel.Dataset.*;

import java.util.Vector;

public class myDataset {

  public static final int REAL = Attribute.REAL;
  public static final int INTEGER = Attribute.INTEGER;
  public static final int NOMINAL = Attribute.NOMINAL;

  private double[][] X = null; //examples array
  private boolean[][] missing = null; //possible missing values
  private int[] outputInteger = null; //output of the data-set as integer values
  private double[] outputReal = null; //output of the data-set as double values
  private String[] output = null; //output of the data-set as string values
  private double[] emax; //max value of an attribute
  private double[] emin; //min value of an attribute
  private double[] weights;
 

  private int nData; // Number of examples
  private int nVars; // Numer of variables
  private int nInputs; // Number of inputs
  private int nClasses; // Number of outputs
  private int[] classList;
  private int[] tipos;
  private String[] clases;
  private String[] variables;
  private String[][] nominales;
  private String[][] nominales_head;

  public InstanceSet IS; //The whole instance set
  private Attribute[] inputs_att;
  private Attribute output_att;

  private double stdev[], average[]; //standard deviation and average of each attribute
  private int instancesCl[];
  private double translation[][];

  public myDataset(myDataset copia, int clase_1, int clase_2) {
    inputs_att = copia.inputs_att.clone();
    output_att = copia.output_att;
    nVars = copia.getnVars();
    nInputs = copia.getnInputs();
    nClasses = 2; //copia.getnClasses();
    clases = copia.clases.clone();
    variables = copia.variables.clone();
    tipos = copia.tipos.clone();
    nominales = new String[nInputs][];
    translation = new double[nInputs][];
    for (int i = 0; i < nominales.length; i++) {
      nominales[i] = copia.nominales[i].clone();
    }
    double[][] X_aux = new double[copia.size()][copia.getnInputs()];
    int[] outputInteger_aux = new int[copia.size()];
    String[] output_aux = new String[copia.size()];
    nData = 0;
    emax = new double[copia.getnInputs()];
    emin = new double[copia.getnInputs()];
    for (int i = 0; i < emax.length; i++) {
      emax[i] = Double.MIN_VALUE;
      emin[i] = Double.MAX_VALUE;
    }
    for (int i = 0; i < copia.size(); i++) {
      if ( (copia.getOutputAsInteger(i) == clase_1) ||
          (copia.getOutputAsInteger(i) == clase_2)) {
        //X_aux[nData] = copia.getExample(i).clone();
        double[] auxiliar = copia.getExample(i).clone();
        for (int j = 0; j < emax.length; j++) {
          X_aux[nData][j] = auxiliar[j];
          if (emax[j] < auxiliar[j]) {
            emax[j] = auxiliar[j];
          }
          if (emin[j] > auxiliar[j]) {
            emin[j] = auxiliar[j];
          }
        }
        outputInteger_aux[nData] = copia.getOutputAsInteger(i);
        output_aux[nData] = copia.getOutputAsString(i);
        nData++;
      }
    }
    nominales_head = new String[nominales.length][];
    for (int i = 0; i < nominales.length; i++) {
      if (tipos[i] == this.NOMINAL) {
        boolean[] auxi = new boolean[nominales[i].length];
        for (int j = 0; j < auxi.length; j++) {
          auxi[j] = false;
        }
        for (int j = 0; j < nData; j++) {
          auxi[ (int) X_aux[j][i]] = true;
        }
        int contador = 0;
        for (int j = 0; j < auxi.length; j++) {
          if (auxi[j]) {
            contador++;
          }
        }
        nominales_head[i] = new String[contador];
        translation[i] = new double[contador + 1];
        contador = 0;
        for (int j = 0; j < auxi.length; j++) {
          if (auxi[j]) {
            nominales_head[i][contador] = nominales[i][j];
            translation[i][contador++] = j;
          }
        }
        translation[i][contador] = -1;
      }
      else {
        translation[i] = new double[1];
        translation[i][0] = 0.0;
        nominales_head[i] = new String[1];
        nominales_head[i][0] = "?";
      }
    }

    X = new double[nData][nInputs];
    outputInteger = new int[nData];
    output = new String[nData];
    for (int i = 0; i < nData; i++) {
      X[i] = X_aux[i].clone();
      outputInteger[i] = outputInteger_aux[i];
      output[i] = output_aux[i];
    }
    classList = new int[2];
    classList[0] = clase_1;
    classList[1] = clase_2;
    copia.computeInstancesPerClass(); //para el cost_sensitive learning

    instancesCl = new int[copia.getnClasses()];
    for (int i = 0; i < instancesCl.length; i++) {
      instancesCl[i] = copia.numberInstances(i);
    }

    weights = new double[copia.getnClasses()];
    if (this.numberInstances(clase_1) < this.numberInstances(clase_2)) {
      weights[clase_1] = 1.0;
      weights[clase_2] = (double)this.numberInstances(clase_1) /
          this.numberInstances(clase_2);
    }
    else {
      weights[clase_1] = (double)this.numberInstances(clase_2) /
          this.numberInstances(clase_1);
      weights[clase_2] = 1.0;
    }
  }

  public myDataset(myDataset copia, int positiva) {
	  inputs_att = copia.inputs_att.clone();
	  output_att = copia.output_att;
	  nVars = copia.getnVars();
	  nInputs = copia.getnInputs();
	  nClasses = 2; //copia.getnClasses();
	  clases = copia.clases.clone();
	  variables = copia.variables.clone();
	  tipos = copia.tipos.clone();
	  nominales = new String[nInputs][];
	  translation = new double[nInputs][];
	  for (int i = 0; i < nominales.length; i++) {
		  nominales[i] = copia.nominales[i].clone();
	  }
	  emax = copia.getemax().clone();
	  emin = copia.getemin().clone();
	  for (int i = 0; i < nominales.length; i++) {
		  //nominales[i] = copia.nominales[i].clone();
		  translation[i] = new double[nominales[i].length+1];
		  for (int j = 0; j < translation[i].length-1; j++){
			  translation[i][j] = 1.0 * j;
		  }
		  translation[i][translation[i].length-1] = -1;
	  }
	  nData = copia.size();
	  nominales_head = nominales.clone();
	  X = new double[nData][nInputs];
	  X = copia.getX().clone();
	  outputInteger = new int[nData];
	  output = new String[nData];
	  int positivos = 0;
	  for (int i = 0; i < nData; i++){
		  outputInteger[i] = 1;
		  output[i] = "negative";
		  if (copia.getOutputAsInteger(i) == positiva){
			  positivos++;
			  outputInteger[i] = 0;
			  output[i] = "positive";
		  }
	  }
	  classList = new int[2];
	  classList[0] = 0;//clase_1;
	  classList[1] = 1;//clase_2;
	  instancesCl = new int[2];
	  instancesCl[0] = positivos;
	  instancesCl[1] = nData - positivos;
	  if (instancesCl[1] < 0) instancesCl[1] = positivos - nData;
	  clases[0] = "positive";
	  clases[1] = "negative";
	  //this.computeInstancesPerClass();
	  weights = new double[2];
	  if (this.numberInstances(0) < this.numberInstances(1)){
		  weights[0] = 1.0;
		  weights[1] = (double) this.numberInstances(0) / this.numberInstances(1);
	  }else{
		  weights[0] = (double) this.numberInstances(1) / this.numberInstances(0);
		  weights[1] = 1.0;
	  }  
  }
  
  
  /**
   * It sets the corresponding classes of the data-set in the global list of classes
   * @param nClasses int total number of classes in the problem
   * @param clases String [] an array with the class labels
   * @param variables String [] an array with the name of the variables
   * @param clase_1 int the first class
   * @param clase_2 int the second class
   */
  void assignClasses(int nClasses, String[] clases, String[] variables,
                     int clase_1, int clase_2) {
    this.nClasses = nClasses;
    classList = new int[2];
    classList[0] = clase_1;
    classList[1] = clase_2;
    for (int i = 0; i < nData; i++) {
      if (outputInteger[i] == 0) {
        outputInteger[i] = clase_1;
      }
      else {
        outputInteger[i] = clase_2;
      }
    }
    instancesCl = new int[nClasses];
    for (int i = 0; i < this.getnData(); i++) {
      instancesCl[this.outputInteger[i]]++;
    }
    this.clases = clases.clone();
    this.variables = variables.clone();
    //System.err.println("Mira -> "+variables.length);
  }

  public String [] getVariablesNames(){
    return this.variables;
  }

  /**
   * It changes the X values for the examples to the real nominal values
   */
  void translateDataSet() {
    for (int i = 0; i < nData; i++) {
      for (int j = 0; j < nInputs; j++) {
        if (this.getTipo(j) == this.NOMINAL) {
          X[i][j] = translation[j][ (int) X[i][j]];
        }
      }
    }
  }

  /**
   * It returns the real nominal value for a given value
   * @param i int Variable
   * @param j int "Position", the nominal value
   * @return int the translation
   */
  public int translation(int i, int j){
    return (int)translation[i][j];
  }

  public int numNominales(int var){
    return translation[var].length;
  }

  double getPeso(int clase) {
    return weights[clase];
  }

  /**
   * Init a new set of instances
   */
  public myDataset() {
    IS = new InstanceSet();
  }

  /**
   * Outputs an array of examples with their corresponding attribute values.
   * @return double[][] an array of examples with their corresponding attribute values
   */
  public double[][] getX() {
    return X;
  }

  /**
   * Output a specific example
   * @param pos int position (id) of the example in the data-set
   * @return double[] the attributes of the given example
   */
  public double[] getExample(int pos) {
    return X[pos];
  }

  /**
   * Returns the output of the data-set as integer values
   * @return int[] an array of integer values corresponding to the output values of the dataset
   */
  public int[] getOutputAsInteger() {
    int[] output = new int[outputInteger.length];
    for (int i = 0; i < outputInteger.length; i++) {
      output[i] = outputInteger[i];
    }
    return output;
  }

  /**
   * Returns the output of the data-set as real values
   * @return double[] an array of real values corresponding to the output values of the dataset
   */
  public double[] getOutputAsReal() {
    double[] output = new double[outputReal.length];
    for (int i = 0; i < outputReal.length; i++) {
      output[i] = outputInteger[i];
    }
    return output;
  }

  /**
   * Returns the output of the data-set as nominal values
   * @return String[] an array of nomianl values corresponding to the output values of the dataset
   */
  public String[] getOutputAsString() {
    String[] output = new String[this.output.length];
    for (int i = 0; i < this.output.length; i++) {
      output[i] = this.output[i];
    }
    return output;
  }

  /**
   * It returns the output value of the example "pos"
   * @param pos int the position (id) of the example
   * @return String a string containing the output value
   */
  public String getOutputAsString(int pos) {
    return output[pos];
  }

  /**
   * It returns the output value of the example "pos"
   * @param pos int the position (id) of the example
   * @return int an integer containing the output value
   */
  public int getOutputAsInteger(int pos) {
    return outputInteger[pos];
  }

  /**
   * It returns the output value of the example "pos"
   * @param pos int the position (id) of the example
   * @return double a real containing the output value
   */
  public double getOutputAsReal(int pos) {
    return outputReal[pos];
  }

  /**
   * Devuelve el valor nominal correspondiente a la clase con valor numerico "clase"
   * @param clase int
   * @return String
   */
  public String nombreClase(int clase) {
    return clases[clase]; //this.output_att.getNominalValue(clase); //Attributes.getOutputAttribute(0).getNominalValue(clase);
  }  
  
  /**
   * It returns an array with the maximum values of the attributes
   * @return double[] an array with the maximum values of the attributes
   */
  public double[] getemax() {
    return emax;
  }

  /**
   * It returns an array with the minimum values of the attributes
   * @return double[] an array with the minimum values of the attributes
   */
  public double[] getemin() {
    return emin;
  }

  public double getMax(int variable) {
    return emax[variable];
  }

  public double getMin(int variable) {
    return emin[variable];
  }

  /**
   * It gets the size of the data-set
   * @return int the number of examples in the data-set
   */
  public int getnData() {
    return nData;
  }

  /**
   * It gets the number of variables of the data-set (including the output)
   * @return int the number of variables of the data-set (including the output)
   */
  public int getnVars() {
    return nVars;
  }

  /**
   * It gets the number of input attributes of the data-set
   * @return int the number of input attributes of the data-set
   */
  public int getnInputs() {
    return nInputs;
  }

  /**
   * It gets the number of output attributes of the data-set (for example number of classes in classification)
   * @return int the number of different output values of the data-set
   */
  public int getnClasses() {
    return nClasses;
  }

  /**
   * This function checks if the attribute value is missing
   * @param i int Example id
   * @param j int Variable id
   * @return boolean True is the value is missing, else it returns false
   */
  public boolean isMissing(int i, int j) {
    return missing[i][j];
  }

  /**
   * It reads the whole input data-set and it stores each example and its associated output value in
   * local arrays to ease their use.
   * @param datasetFile String name of the file containing the dataset
   * @param train boolean It must have the value "true" if we are reading the training data-set
   * @throws IOException If there ocurs any problem with the reading of the data-set
   */
  public void readClassificationSet(String datasetFile, boolean train) throws
      IOException {
    try {
      // Load in memory a dataset that contains a classification problem
      IS = new InstanceSet();
      IS.readSet(datasetFile, train);
      IS.setAttributesAsNonStatic();
      inputs_att = IS.getAttributeDefinitions().getInputAttributes();
      output_att = IS.getAttributeDefinitions().getOutputAttribute(0);
      
      nData = IS.getNumInstances();
      nInputs = IS.getAttributeDefinitions().getInputNumAttributes();
      nVars = nInputs + IS.getAttributeDefinitions().getOutputNumAttributes();

      // outputIntegerheck that there is only one output variable
      if (IS.getAttributeDefinitions().getOutputNumAttributes() > 1) {
        System.out.println(
            "This algorithm can not process MIMO datasets");
        System.out.println(
            "All outputs but the first one will be removed");
        System.exit(1);
      }
      boolean noOutputs = false;
      if (IS.getAttributeDefinitions().getOutputNumAttributes() < 1) {
        System.out.println(
            "This algorithm can not process datasets without outputs");
        System.out.println("Zero-valued output generated");
        noOutputs = true;
        System.exit(1);
      }

      // Initialice and fill our own tables
      X = new double[nData][nInputs];
      missing = new boolean[nData][nInputs];
      outputInteger = new int[nData];
      outputReal = new double[nData];
      output = new String[nData];

      // Maximum and minimum of inputs
      emax = new double[nInputs];
      emin = new double[nInputs];
      for (int i = 0; i < nInputs; i++) {
        if (inputs_att[i].getType() == this.NOMINAL){
          emin[i] = 0.0;
          emax[i] = (double)inputs_att[i].getNumNominalValues();
        }else{
          emax[i] = inputs_att[i].getMaxAttribute();
          emin[i] = inputs_att[i].getMinAttribute();
        }
        //System.out.println("Att["+i+"]: ("+emin[i]+","+emax[i]+")");
      }
      // All values are casted into double/integer
      nClasses = 0;
      for (int i = 0; i < nData; i++) {
        Instance inst = IS.getInstance(i);
        for (int j = 0; j < nInputs; j++) {
          X[i][j] = IS.getInputNumericValue(i, j); //inst.getInputRealValues(j);
          missing[i][j] = inst.getInputMissingValues(j);
          if (missing[i][j]) {
            X[i][j] = emin[j] - 1;
          }
        }

        if (noOutputs) {
          outputInteger[i] = 0;
          output[i] = "";
        }
        else {
          outputInteger[i] = (int) IS.getOutputNumericValue(i, 0);
          output[i] = IS.getOutputNominalValue(i, 0);
        }
        if (outputInteger[i] > nClasses) {
          nClasses = outputInteger[i];
        }
      }
      nClasses++;
      System.out.println("Number of classes=" + nClasses);

    }
    catch (Exception e) {
      System.out.println("DBG: Exception in readSet");
      e.printStackTrace();
    }
    classList = new int[nClasses];
    for (int i = 0; i < nClasses; i++) {
      classList[i] = i;
    }
    computeStatistics();
    computeInstancesPerClass();
    computeWeightsProportional();
    variables = new String[nVars];
    clases = new String[nClasses];
    tipos = new int[nInputs];
    nominales = new String[nInputs][];
    nominales_head = new String[nInputs][];
    for (int i = 0; i < nInputs; i++) {
      variables[i] = inputs_att[i].getName();
      tipos[i] = inputs_att[i].getType();
      if (inputs_att[i].getNumNominalValues() > 0) {
        nominales[i] = new String[inputs_att[i].getNumNominalValues()];
        nominales_head[i] = new String[inputs_att[i].getNumNominalValues()];
        for (int j = 0; j < nominales[i].length; j++) {
          nominales[i][j] = inputs_att[i].getNominalValue(j);
          nominales_head[i][j] = inputs_att[i].getNominalValue(j);
        }
      }
      else {
        nominales[i] = new String[1];
        nominales[i][0] = "?";
        nominales_head[i] = new String[1];
        nominales_head[i][0] = "?";
      }
    }
    variables[nInputs] = output_att.getName();
    for (int i = 0; i < nClasses; i++) {
      clases[i] = output_att.getNominalValue(i);
    }
  }

  /**
   * It transform the input space into the [0,1] range
   */
  public void normalize() {
    int atts = this.getnInputs();
    double maxs[] = new double[atts];
    for (int j = 0; j < atts; j++) {
      maxs[j] = 1.0 / (emax[j] - emin[j]);
    }
    for (int i = 0; i < this.getnData(); i++) {
      for (int j = 0; j < atts; j++) {
        if (isMissing(i, j)) {
          ; //this process ignores missing values
        }
        else {
          X[i][j] = (X[i][j] - emin[j]) * maxs[j];
        }
      }
    }
  }

  /**
   * It checks if the data-set has any missing value
   * @return boolean True if it has some missing values, else false.
   */
  public boolean hasMissingAttributes() {
    return (this.sizeWithoutMissing() < this.getnData());
  }

  /**
   * It return the size of the data-set without having account the missing values
   * @return int the size of the data-set without having account the missing values
   */
  public int sizeWithoutMissing() {
    int tam = 0;
    for (int i = 0; i < nData; i++) {
      int j;
      for (j = 1; (j < nInputs) && (!isMissing(i, j)); j++) {
        ;
      }
      if (j == nInputs) {
        tam++;
      }
    }
    return tam;
  }

  public int size() {
    return nData;
  }

  /**
   * It computes the average and standard deviation of the input attributes
   */
  private void computeStatistics() {
    stdev = new double[this.getnVars()];
    average = new double[this.getnVars()];

    for (int i = 0; i < this.getnInputs(); i++) {
      average[i] = 0;
      for (int j = 0; j < this.getnData(); j++) {
        if (!this.isMissing(j, i)) {
          average[i] += X[j][i];
        }
      }
      average[i] /= this.getnData();
    }
    average[average.length - 1] = 0;
    for (int j = 0; j < outputReal.length; j++) {
      average[average.length - 1] += outputReal[j];
    }
    average[average.length - 1] /= outputReal.length;

    for (int i = 0; i < this.getnInputs(); i++) {
      double sum = 0;
      for (int j = 0; j < this.getnData(); j++) {
        if (!this.isMissing(j, i)) {
          sum += (X[j][i] - average[i]) * (X[j][i] - average[i]);
        }
      }
      sum /= this.getnData();
      stdev[i] = Math.sqrt(sum);
    }

    double sum = 0;
    for (int j = 0; j < outputReal.length; j++) {
      sum += (outputReal[j] - average[average.length - 1]) *
          (outputReal[j] - average[average.length - 1]);
    }
    sum /= outputReal.length;
    stdev[stdev.length - 1] = Math.sqrt(sum);
  }

  /**
   * It return the standard deviation of an specific attribute
   * @param position int attribute id (position of the attribute)
   * @return double the standard deviation  of the attribute
   */
  public double stdDev(int position) {
    return stdev[position];
  }

  /**
   * It return the average of an specific attribute
   * @param position int attribute id (position of the attribute)
   * @return double the average of the attribute
   */
  public double average(int position) {
    return average[position];
  }

  public void computeInstancesPerClass() {
    instancesCl = new int[nClasses];
    for (int i = 0; i < this.getnData(); i++) {
      instancesCl[this.outputInteger[i]]++;
    }
  }

  public int numberInstances(int clas) {
    return instancesCl[clas];
  }

  public int numberValues(int attribute) {
    return (int) emax[attribute];
    //return inputs_att[attribute].getNumNominalValues();
  }

  public String getOutputValue(int intValue) {
    return clases[intValue];
  }

  public int getOutputValue(String value) {
    for (int i = 1; i < nClasses; i++) {
      if (clases[i].compareToIgnoreCase(value) == 0) {
        return i;
      }
    }
    return 0;
  }

  public int getAttribute(String att) {
    for (int i = 0; i < nVars; i++) {
      if (variables[i].equalsIgnoreCase(att)) {
        return i;
      }
    }
    return -1;
  }

  public int getTipo(int variable) {
    return tipos[variable];
  }

  /**
   * It return the universe of discourse of the variable
   * @return double[][] The minimum and maximum range of each variable
   */
  public double[][] devuelveRangos() {
    double[][] rangos = new double[this.getnVars()][2];
    for (int i = 0; i < this.getnInputs(); i++) {
      rangos[i][0] = emin[i];
      rangos[i][1] = emax[i];
    }
    //rangos[this.getnVars() - 1][0] = output_att.getMinAttribute();
    //rangos[this.getnVars() - 1][1] = output_att.getMaxAttribute();
    rangos[this.getnVars() - 1][0] = 0;
    rangos[this.getnVars() - 1][1] = nClasses - 1;
    return rangos;
  }

  /**
   * It computes an array that stores a boolean value that indicates if a given attribute
   * has nominal type or not
   * @return boolean [] an array that for each position tells if the attribute is nominal
   */
  public boolean[] getNominals() {
    boolean[] nominals = new boolean[this.nInputs];
    for (int i = 0; i < nInputs; i++) {
      nominals[i] = (this.getTipo(i) == this.NOMINAL);
    }
    return nominals;
  }

  /**
   * It returns an array with the names of the input attributes
   * @return String[] an array with the names of the input attributes
   */
  public String[] varNames() {
    String nombres[] = new String[nInputs];
    for (int i = 0; i < nInputs; i++) {
      //nombres[i] = inputs_att[i].getName();
      nombres[i] = variables[i];
    }
    return nombres;
  }

  /**
   * It returns an array with the labels of the output class
   * @return String[] an array with the labels of the output class
   */
  public String[] classNames() {
    String clases_[] = new String[nClasses];
    for (int i = 0; i < nClasses; i++) {
      clases_[i] = clases[i]; //output_att.getNominalValue(i);
    }
    return clases_;
  }

  private void computeWeightsProportional() {
    double max = 0.0;
    weights = new double[getnClasses()];
    for (int i = 0; i < weights.length; i++) {
      if (numberInstances(i) > 0) {
        weights[i] = 1.0 * size() / numberInstances(i);
        if (weights[i] > max) {
          max = weights[i];
        }
      }
      else {
        weights[i] = 0.0;
      }
    }
    for (int i = 0; i < weights.length; i++) {
      weights[i] /= max;
    }

  }

  public boolean empty() {
    if ( (this.numberInstances(classList[0]) == 0) ||
        (this.numberInstances(classList[1]) == 0)) {
      return true;
    }
    return false;
  }

  public double ir(){
	  if (this.numberInstances(classList[0]) >
	  this.numberInstances(classList[1])) {
		  return ( (double)this.numberInstances(classList[0]) /
				  this.numberInstances(classList[1]));
	  }
	  else {
		  return ( (double)this.numberInstances(classList[1]) /
				  this.numberInstances(classList[0]));
	  }
  }
  
  public boolean noBalanceado() {
    //System.err.println(numberInstances(listaClases[0])+" instancias vs "+numberInstances(listaClases[1]));
    return ir() > 1.5;
  }

  public int n_minoritaria() {
    if (this.numberInstances(classList[0]) >
        this.numberInstances(classList[1])) {
      return this.numberInstances(classList[1]);
    }
    else {
      return this.numberInstances(classList[0]);
    }

  }

  public String doHeader() {
    String cadena = new String("");
    cadena += "@relation unknown\n";
    for (int i = 0; i < this.nInputs; i++) {
      //Attribute a = inputs_att[i];
      cadena += "@attribute " + variables[i];
      if (tipos[i] == this.INTEGER) {
        cadena += " integer [" + (int) emin[i] + "," + (int) emax[i] + "]\n";
      }
      else if (tipos[i] == this.REAL) {
        cadena += " real [" + emin[i] + "," + emax[i] + "]\n";
      }
      else {
        cadena += " {";
        int j;
        for (j = 0; j < nominales_head[i].length - 1; j++) {
          cadena += nominales_head[i][j] + ",";
        }
        cadena += nominales_head[i][j] + "}\n";
      }
    }
    //Attribute a = output_att;
    cadena += "@attribute " + variables[nInputs];
    cadena += " {";
    int i = 0;
    for (; i < nClasses-1; i++){
    	cadena += clases[classList[i]] + ", ";
    }
    cadena += clases[classList[i]] + "}\n";

    cadena += "@data\n";
    return cadena;
  }

  public String printDataSet() {
    String cadena = new String("");
    /*if (!preprocesado){
      cadena += doHeader();
    }else{
      cadena += this.copyHeader();
    }*/
    cadena += doHeader();
    //Attributes.clearAll();
    for (int i = 0; i < size(); i++) {
      double[] ejemplo = this.getExample(i);
      for (int j = 0; j < this.getnInputs(); j++) {
        if (this.getTipo(j) == this.NOMINAL) {
          cadena += nominales[j][ (int) ejemplo[j]] + ", ";
        }
        else if (this.getTipo(j) == this.INTEGER) {
          cadena += (int) ejemplo[j] + ", ";
        }
        else {
          cadena += ejemplo[j] + ", ";
        }
      }
      cadena += this.getOutputAsString(i) + "\n";
    }
    return cadena;
  }

  public String claseMasFrecuente(){
    /*
	int [] clases = new int[nClasses];
    for (int i = 0; i < this.outputInteger.length; i++){
      clases[outputInteger[i]]++;
    }*/
    int claseMayoritaria = 0;
    for (int i = 1; i < nClasses; i++){
      if (instancesCl[claseMayoritaria] < instancesCl[i]){
        claseMayoritaria = i;
      }
    }
    return this.getOutputValue(claseMayoritaria);
  }

  public double valorReal(int atributo, String valorNominal){
    Vector nominales = inputs_att[atributo].getNominalValuesList();
    int aux = nominales.indexOf(valorNominal);
    return 1.0*aux;
  }

  public int claseNumerica(String valorNominal){
    Vector nominales = output_att.getNominalValuesList();
    int aux = nominales.indexOf(valorNominal);
    if ((nominales.size() == 2)&&(aux >= 0)){
      aux = classList[aux];
    }
    return aux;
  }

  public String valorNominal(int atributo, double valorReal){
    Vector nominales = inputs_att[atributo].getNominalValuesList();
    return (String)nominales.get((int)valorReal);
  }

  public int totalNominales(int atributo){
    return inputs_att[atributo].getNumNominalValues();
  }

}
