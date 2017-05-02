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

public class RuleBase {

  ArrayList<Rule> ruleBase;
  myDataset train;

  public RuleBase() {
    ruleBase = new ArrayList<Rule> ();
  }

  /**
   * Obtengo la base de reglas a traves del fichero de reglas (extraido a partir del arbol de decision)
   * @param reglas String
   * @param train myDataset conjunto de datos de entrenamiento
   */
  public RuleBase(myDataset train, String reglas) {
    ruleBase = new ArrayList<Rule> ();
    this.train = train;
    StringTokenizer tokens = new StringTokenizer(reglas, "\n");
    while (tokens.hasMoreTokens()) {
      String regla = tokens.nextToken();
      //System.err.println("Regla -> "+regla);
      Rule r = new Rule(train, regla);
      ruleBase.add(r);
    }
  }

  /**
   * Constructor
   * @param regla Default rule
   * @param train myDataset training data
   */
  public RuleBase(myDataset train, Rule regla) {
    ruleBase = new ArrayList<Rule> ();
    this.train = train;
    ruleBase.add(regla);
  }  
  
  public String printString() {
    String cadena = new String("");
    cadena += "Number of Rules: " + ruleBase.size() + "\n";
    for (int i = 0; i < ruleBase.size(); i++) {
      cadena += "Rule[" + (i + 1) + "]: " + ruleBase.get(i).printString();
    }
    return cadena;
  }

  public int size() {
    return ruleBase.size();
  }

  /**
   * Detect those rules that cover an small disjunct
   */
  public void coverExamples() {
    for (int i = 0; i < this.size(); i++) {
      ruleBase.get(i).coverExamples();
    }
  }

}
