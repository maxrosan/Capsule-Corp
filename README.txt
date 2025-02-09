       ==================================================================
			   Universidade de S�o Paulo
		      Instituto de Matem�tica e Estat�stica
       ==================================================================
        Nome: Max Rosan dos Santos Junior             NUSP: 
              Patr�cia Ara�jo de Oliveira             NUSP: 7807127

        MAC5742 - Introdu��o � Computa��o Paralela e Distribu�da
        Professor: Marco Dimas Gubitoso
       ==================================================================

		        --------------------------
                       Projeto de C�psula Alien�gena
		        --------------------------

-----------
Objetivo
-----------
O programa tem como objetivo simular a distribui��o de calor em pastilhas de uma
c�psula que est� entrando na atmosfera em alta velocidade. Seguindo especifica-
��o do problema, s�o consideradas diversas simplifica��es e certas complica��es.
Como forma de aplica��o dos conhecimentos adquiridos, a execu��o do programa
deve ocorrer de forma paralelizada, aproveitando recursos disponibilizados pelos
processadores multi-core.

----------------------
Considera��es Gerais 
----------------------

A fim de n�o alterar posi��o de todas as pastilhas em rela��o ao vetor posi��o, 
rotacionamos os vetores posi��o e velocidade de modo "s�ncrono",fazendo com que o
vetor posi��o aponte para baixo (<0, 0, z>, z<0) e mantendo a posi��o relativa
entre esses vetores, n�o alterando, assim, resultado dos c�lculos.

Cada passo da simula��o foi dividido em duas partes: c�lculo da temperatura e 
atualiza��o da temperatura. O c�lculo da temperatura de cada pastilha � 
realizado com base nos valores obtidos no passo anterior, o novo valor encontrado
n�o ser� utilizado nos c�lculos da temperatura das pastilhas vizinhas. Na
segunda parte, a atualiza��o da temperatura, o valor encontrado na parte descrita
anteriormente � definido como valor corrente da temperatura para ent�o ser utilizado
no pr�ximo passo. Ent�o cada passo da simula��o pode ser dividido em dois la�os: o
loop para c�lculo da temperatura e o para atualiza��o dos valores. Como o c�lculo
da temperatura da pastilha s� usa valores do passo anterior, � poss�vel paralelizar
as itera��es; assim, al�m da paraleliza��o realizada nesses la�os, tamb�m foi realizado
unrolling.

----------------------
Entradas do Programa
----------------------
O programa recebe como entrada os par�metros que descrevem a c�psula, tanto geo-
metria quanto propriedades das pastilhas que a revestem, e o ambiente, na se-
guinte ordem:

1. h            - corresponde � altura da c�psula, medida da base at� o v�rtice.
2. a            - corresponde ao fator de forma o paraboloide.
3. d            - corresponde ao lado da pastilha.
4. alpha        - corresponde a um par�metro da fun��o de atrito.
5. t_0          - corresponde a um par�metro da fun��o de atrito.
6. delta        - corresponde ao par�metro da fun��o de dissipa��o.
7. theta_crit   - corresponde � temperatura na qual a pastilha se desintegra.
8. theta_0      - corresponde � temperatura inicial.
9. pos          - Vetor posi��o (valores de x, y e z separados por espa�o).
10.vel          - Vetor velocidade (valores de x, y e z separados por espa�o).
11.steps        - corresponde ao n�mero de itera��es que devem ser calculadas.

O programa faz a leitura dessas vari�veis atrav�s do arquivo 'entrada.txt', o qual
possui, nessa ordem, apenas os valores que correspondem a cada uma delas.

----------------
Como executar?
----------------
Para facilitar esse processo foi criado o 'Makefile', para compilar o programa.
Para isto basta entrar no diret�rio,onde se encontram os arquivos do programa, e
digitar o comando 'make', utilizando o terminal.

Ao compilar o c�digo fonte, ser� criado um arquivo execut�vel de nome 'ep'.
� importante verificar se o arquivo de entrada 'entrada.txt' est� presente no
diret�rio. Para modificar os valores que deseja passar como entrada, basta
editar este arquivo.

Para executar, execute no terminal "./ep".

---------
Sa�das
---------
O programa simula a distribui��o das temperaturas a cada passo e escreve um
arquivo texto na sa�da (saida.txt). Al�m disso, o programa imprime no terminal a
temperatura m�dia das pastilhas e do rejunte.

-----------------------------
Formato do arquivo de sa�da
-----------------------------
O arquivo de sa�da (saida.txt) possui os valores colocados em linhas, da
seguinte forma:

-> Primeira linha: a, h e d.
-> Segunda linha: temperatura da calota.
-> Terceira linha: n�mero de an�is.
-> Linhas subsequentes: um anel por linha, com o n�mero de pastilhas na primeira
posi��o, seguido pelas temperaturas de cada uma. Come�ando pelo topo do anel.

OBS: Temperatura positiva indica a temperatura da pastilha existente na posi��o,
enquanto uma negativa indica a temperatura do rejunte na posi��o.
Os valores das temperaturas s�o dados em Celsius.


               ==================================================
                          Data: 28 de Novembro de 2011
               ==================================================

