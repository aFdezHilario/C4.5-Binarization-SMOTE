package multi_C45;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import java.util.ArrayList;
import java.util.StringTokenizer;

public class Rule {

  ArrayList<Selector> antecedent;
  String clase;
  myDataset train;
  int coveredExamples[];
  int examplesCoveredOK[];
  int nCovered,nCoveredOk; //number of covered examples
  double fitness; //rule obtained by GA
  int codigoRegla; //rule produced by GA

  public Rule() {
    antecedent = new ArrayList<Selector> ();
    coveredExamples = new int[1];
  }

  public Rule(String clase, myDataset train) {
    antecedent = new ArrayList<Selector> ();
    this.train = train;
    this.clase = clase;
    coveredExamples = new int[train.size()];
    examplesCoveredOK = new int[train.size()];
  }

  public Rule(myDataset train, String linea) {
    antecedent = new ArrayList<Selector> ();
    this.train = train;
    coveredExamples = new int[train.size()];
    examplesCoveredOK = new int[train.size()];
    String [] nombres = train.varNames();
    StringTokenizer campo = new StringTokenizer(linea, " ");
    campo.nextToken(); //RULE-X:
    String aux = campo.nextToken(); //IF
    while(!aux.equalsIgnoreCase("THEN")){
      String atributo = campo.nextToken();
      String operador = campo.nextToken();
      String valor = campo.nextToken();
      Selector s = new Selector(atributo,operador,valor,train);
      s.adjuntaNombres(nombres);
      antecedent.add(s);
      aux = campo.nextToken();
    }
    campo.nextToken(); //class
    campo.nextToken(); //=
    clase = campo.nextToken();
  }

  public void addSelector(Selector s) {
    antecedent.add(s);
  }

  public String printString() {
    String cadena = new String("");
    cadena += "IF ";
    if (antecedent.size()>0){
	    for (int i = 0; i < antecedent.size()-1; i++) {
	      cadena += antecedent.get(i).printString()+ "AND ";
	    }
	    cadena += antecedent.get(antecedent.size()-1).printString();
    }else{
    	cadena += " * ";
    }
    cadena += " THEN Class = " + clase + " ("+nCoveredOk+"/"+nCovered+")\n";
    return cadena;
  }

  public Rule copy(){
    Rule r = new Rule(clase, train);
    r.antecedent = new ArrayList<Selector>();
    for (int i = 0; i < antecedent.size(); i++){
      r.antecedent.add(antecedent.get(i).copia());
    }
    r.nCovered = nCovered;
    r.nCoveredOk = nCoveredOk;
    r.coveredExamples = new int[coveredExamples.length];
    r.coveredExamples = coveredExamples.clone();
    r.examplesCoveredOK = new int[examplesCoveredOK.length];
    r.examplesCoveredOK = examplesCoveredOK.clone();
    r.fitness = fitness;
    r.codigoRegla = codigoRegla;
    return r;
  }

  public int covered(){
    return nCovered;
  }

  public int coveredOK(){
    return nCoveredOk;
  }

  public void coverExamples(){
    nCovered = nCoveredOk = 0;
    for (int i = 0; i < train.size(); i++){
      double [] ejemplo = train.getExample(i);
      if (this.cubre(ejemplo)){
        coveredExamples[nCovered] = i;
        nCovered++;
        if (train.getOutputAsString(i).compareToIgnoreCase(this.clase) == 0){
          examplesCoveredOK[nCoveredOk] = i;
          nCoveredOk++;
        }
      }
    }
  }
  
  public double confidence(){
	  return 1.0*nCoveredOk/nCovered;
  }

  public boolean cubre(double [] ejemplo){
    boolean cubierto = true;
    for (int i = 0; (i < antecedent.size())&&(cubierto); i++){
      cubierto = cubierto && (antecedent.get(i).cubre(ejemplo));
    }
    return cubierto;
  }

  public int size(){
    return antecedent.size();
  }

  public boolean contieneAtributo(int att){
    boolean contiene = false;
    for (int i = 0; i < antecedent.size() && !contiene; i++){
      contiene = (antecedent.get(i).atributo == att);
    }
    return contiene;

  }

}
