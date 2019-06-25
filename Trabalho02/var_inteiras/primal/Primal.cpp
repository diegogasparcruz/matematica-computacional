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
            for(int j = 0;j < qtdDestinos; j++)
                entrada >> C[i][j];

        //lendo valores de oferta
        for(int i = 0 ;i < qtdOrigens; i++)
            entrada >> origem[i];

        //lendo valores de demanda
        for(int j = 0; j < qtdDestinos; j++)
            entrada >> destino[j];


        //fazendo verificação pra ver se tudo foi lido corretamente
        /*entrada >> ultimaLinha;
        if(ultimaLinha == 999)
            cout << "Entrada lida com sucesso" << endl;
        else{
            cout << "Erro ao ler entrada" << endl;
            return 1;
        }*/


        //restrição de não negatividade pra X
        IloArray<IloNumVarArray> x(env, qtdOrigens);
        for (int i = 0; i < qtdOrigens; i++){
            x[i] = IloNumVarArray(env, qtdDestinos, 0, IloInfinity, ILOINT);
        }

       
        //definindo função objetiva
        IloExpr Z(env);
        for(int i = 0; i < qtdOrigens; i++)
            for(int j = 0; j < qtdDestinos; j++){
                Z += C[i][j] * x[i][j];
            }


        //definindo tipo da função
        modelo.add(IloMinimize(env, Z));


        //adicionando restrições de fornecimento (valor enviado deve ser menor ou igual ao máximo disponível)
        IloExpr Expr(env);
        for(int i = 0; i < qtdOrigens; i++){
            for(int j = 0; j < qtdDestinos; j++){
                Expr += x[i][j];
            }
            modelo.add(Expr <= origem[i]);
            Expr.end();
            Expr = IloExpr(env);
        }


        //adicionando restrições de demanda (valor recebido deve ser maior ou igual ao que o destino precisa)
        IloExpr Expr2(env);
        for(int j = 0; j < qtdDestinos; j++){
            for(int i = 0;i < qtdOrigens; i++){
                Expr2 = Expr2 + x[i][j];
            }
            modelo.add(Expr2 >= destino[j]);
            Expr2.end();
            Expr2 = IloExpr(env);
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
