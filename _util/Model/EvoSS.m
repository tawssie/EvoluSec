function Model=EvoSS(Algo)
% Create a data structure containing evolutionary model of secondary structure
% 7 labels of "BEGHIST" represent secondary structures
% B = beta bridge
% E = beta strand
% G = 3-10 helix
% H = alpha helix
% I = pi helix
% S = Bend
% T = Turn

AlgoOption={'DSO-SS' 'JTT-SS'};
% default Algo for model is DSO-SS
if ~exist('Algo','var')
    Algo='DSO-SS';
end
    
Algo=cellstr(Algo);

Model.AAOrder='BEGHIST';
Model.F=[]; % normalised frequencies
Model.Q=[]; % transistion (instaneous) rate matrix
Model.V=[]; % eigenvectors
Model.D=[]; % eigenvalues

if seqmatch(Algo,AlgoOption)==2
    [Model.F, Model.Q, Model.V, Model.D]=GetFRVD4JTT_SS;
else
    [Model.F, Model.Q, Model.V, Model.D]=GetFRVD4DSO_SS;
end









function [F, Q, V, D]=GetFRVD4DSO_SS

% normalised frequencies
F=[0.01014266893169180
   0.26401205821245700
   0.05599759650153780
   0.42824733005847200
   0.00005733475407192
   0.10028203535081600
   0.14126097619095400];


% 1 PAM for SS model based on DSO approach
OnePAM=[0.976623887403790       0.0176994108098560      0                       0                       0                       0.00470992724632022     0.000966774540034149
        0.000626636369881649	0.998063044162005       6.67007900714360e-05	2.10634073909798e-05	0                       0.000986469579477554	0.000236085691173899
        0                       0.000405244272098312	0.972032813063742       0.00777429143051762     0                       0.00331127227596121     0.0164763789576814
        0                       1.62520241877680e-05	0.000987310469406908	0.994923273944346       4.87560725633041e-05	0.000538348301219816	0.00348605918827624
        0                       0                       0                       0.213467577664294       0.644220703892843       0                       0.142311718442863
        0.000520346400879469	0.00307825976099223     0.00170071113129553     0.00217723888789041     0                       0.980150154139082       0.0123732896798602      
        7.37030248644819e-05	0.000508361889449888	0.00583954735464741     0.00972879928211161     9.07114152178238e-05	0.00853821195737767     0.975220665076331];


[V_1pam,D_1pam]=eig(OnePAM);

D=diag(log(diag(D_1pam))*100);
% Q = transition rate matrix (also called instantaneous rate matrix)
Q=V_1pam*D/V_1pam; 


% correction to avoid numerical error
Q(eye(size(Q))==1)=0;
Q(Q<0)=0;

Q(eye(size(Q))==1)=-sum(Q,2);

% calculate V and D according to corrected Q
[V,D]=eig(Q);




function [F, Q, V, D]=GetFRVD4JTT_SS


F=[0.0124675060282499
   0.303201526484448
   0.0459107290845415
   0.400735777795490
   0.000182455992094854
   0.103276813588554
   0.134225191026622];

OnePAM=[0.973720660114133    ,0.0189576696231223    ,0.000177648845196594 ,0.000596392551731424 ,2.53784064566564e-05   ,0.00497416766550464  ,0.00154808279385604
        0.000839027762549938 ,0.998006326267034     ,3.31342958436723e-05 ,4.88590464135506e-05 ,0                      ,0.000936745855377039 ,0.000135906772782520
        4.78004119697164e-05 ,0.000201444593300948  ,0.971186594527683    ,0.00759685118804421  ,6.82863028138805e-06   ,0.00380013275159245  ,0.0171603478971282
        1.97554980051595e-05 ,3.65686877967845e-05  ,0.000935233682159144 ,0.994693757301763    ,2.56401144322282e-05   ,0.000701110014310765 ,0.00358793470153279
        0.00131621218726312  ,0                     ,0.00131621218726312  ,0.0401444717115251   ,0.919052950483318      ,0.00197431828089468  ,0.0361958351497358
        0.000616067240671723 ,0.00262142897306233   ,0.00174919091547864  ,0.00262142897306233  ,4.71480031126318e-06   ,0.979724787061465    ,0.0126623820359492
        0.000144760047287844 ,0.000287146979046378  ,0.00596363932515329  ,0.0101284570790904   ,6.52606770559950e-05   ,0.00956009590982094  ,0.973850639982545];


[V_1pam,D_1pam]=eig(OnePAM);

D=diag(log(diag(D_1pam))*100);
% Q = transition rate matrix (also called instantaneous rate matrix)
Q=V_1pam*D/V_1pam;

% correction to avoid numerical error
Q(eye(size(Q))==1)=0;
Q(Q<0)=0;

Q(eye(size(Q))==1)=-sum(Q,2);

% calculate V and D according to corrected Q
[V,D]=eig(Q);











