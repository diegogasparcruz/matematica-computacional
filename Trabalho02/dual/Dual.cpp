#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

int main(int argc, char* argv[]) {

    if(argc < 2){
        cout << "é necessário passar o caminho para um arquivo" << endl;
        return 1;
    }

    IloEnv env;
    
    try {
        IloModel modelo(env);

        IloInt qtdOrigens, qtdDestinos, ultimaLinha;

        //carregando arquivo
        ifstream entrada(argv[1]);
        
        //lendo as dimenções da matriz de coeficientes
        entrada  >> qtdOrigens;
        entrada  >> qtdDestinos;

        //definindo a matriz de custos
        IloArray<IloNumArray> C(env, qtdOrigens);
        for (int i = 0; i < qtdOrigens; i++){
            C[i] = IloNumArray (env, qtdDestinos);
        }


        //definindo array de oferta
        IloNumArray origem(env, qtdOrigens);

        //definindo array de demanda
        IloNumArray destino(env, qtdDestinos);


        //lendo a matriz de custos
        for(int i = 0;i < qtdOrigens; i++)
            for(int j = 0; j < qtdDestinos; j++)
                entrada >> C[i][j];


        //lendo valores de oferta
        for(int i = 0 ;i < qtdOrigens; i++)
            entrada >> origem[i];

        //lendo valores de demanda
        for(int j = 0; j < qtdDestinos; j++)
            entrada >> destino[j];


        IloNumVarArray Y1(env);
        IloNumVarArray Y2(env);

        for(int i = 0; i < qtdOrigens; i++)
            Y1.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
            
        for(int j = 0; j < qtdDestinos; j++)
            Y2.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT)); 
            

        //definindo função objetivo
        IloExpr W(env);
        for(int i = 0; i < qtdOrigens; i++)
            W -= Y1[i] * origem[i];  

        for(int j = 0; j < qtdDestinos; j++)
            W += Y2[j] * destino[j];

        
        //definindo tipo da função
        modelo.add(IloMaximize(env, W));
        

        //adicionando restrições
        for(int i = 0; i < qtdOrigens; i++){
            for (int j = 0; j < qtdDestinos; ++j) {
                IloExpr exp2(env);
                exp2 = -Y1[i] + Y2[j];
                
                modelo.add(exp2 <= C[i][j]);
                exp2.end();
            }
        }
        

        //resolvendo problema
        clock_t start, end;
        double elapsed = 0;
        start = clock();

        IloCplex cplex (env);
        cplex.setOut(env.getNullStream());
        cplex.extract(modelo);
        cplex.exportModel ( "modelo.lp" );

        if ( !cplex.solve() ) {
            env.error() << "Não tem sulução" << endl;
            throw (-1);
        }
        cout << "Problema resolvido!" << endl;
        cout << "Valor da função objetivo = " << cplex.getObjValue() << endl;

        end = clock();
        elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
        cout << "Tempo de execução " << elapsed << "s" << endl;

    } catch (IloException& ex) {
        cerr << "Error Cplex: " << ex << endl;
    } catch (...) {
        cerr << "Error Cpp" << endl;
    }
    env.end();
    return 0;
}
