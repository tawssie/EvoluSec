function Data=Evolusec_BoundaryTest1()
% The script calculates the secondary structure variation around the
% boundaries between modern and actual, reconstructed ancestral protein
% structures.
% The ancestral secondary structure with MSSAs are pre-calculated
% by Evolu-sec-ASSR and embedded here.
% 
% Thus,  call the Data=Evolusec_BoundaryTest1();
%
% Data.Name            % describe the compairson
% Data.TreePairTable   % specific pairwise comparison in the tree
% Data.MSSATable       % multiple secondary structure alignment
% Data.CountTable      % a table count the variation
%
% Well-formated table is in the supplementary file
% "SummaryTable.pdf"


Data=GetDataset();

Data(1).CountTable=[];

for k=1:numel(Data)
    TreePairTable=Data(k).TreePairTable;
    MSSATable=Data(k).MSSATable;
    CountTable=zeros(size(TreePairTable,1),4);
    for i=1:size(TreePairTable,1)

        AncName   =TreePairTable(i,1);
        ModernName=TreePairTable(i,2);

        AncSS=MSSATable(strcmp(MSSATable(:,1),AncName),2);
        ModernSS=MSSATable(strcmp(MSSATable(:,1),ModernName),2);

        % The boundary is defined by secondary parameter,
        %     i.e., ModernSS (the Modern protein)
        [CountTable(i,1),CountTable(i,2),CountTable(i,3),CountTable(i,4)]=...
            CountVariation(AncSS{:},ModernSS{:});
    end
    Data(k).CountTable=CountTable;
end


function [BC,nBC,nBnC,BnC]=CountVariation(AncSS,ModernSS)
% Here, isBoundary is defined by ModernSS (modern proteins)
% BC   = counts for SS in Boundary     with Change
% nBC  = counts for SS in Non-Boundary with Change
% nBnC = counts for SS in Non-Boundary with No Change
% BnC  = counts for SS in Boundary     with No Change

isBoundary=true(size(ModernSS));
for i=2:numel(ModernSS)-1
    
    if (ModernSS(i-1)==ModernSS(i))&&(ModernSS(i)==ModernSS(i+1))
        isBoundary(i)=false;
    end
end

isGap    =AncSS=='-'|ModernSS=='-';
isUnknown=AncSS=='_'|ModernSS=='_';
isChange =AncSS~=ModernSS;

BC  =numel(AncSS( isBoundary&~isGap&~isUnknown&isChange));
nBC =numel(AncSS(~isBoundary&~isGap&~isUnknown&isChange));
nBnC=numel(AncSS(~isBoundary&~isGap&~isUnknown&~isChange));
BnC =numel(AncSS( isBoundary&~isGap&~isUnknown&~isChange));






function Data=GetDataset()

Data.Name=[];
Data.TreePairTable=[];
Data.MSSATable=[];

Data=repmat(Data,12,1);

% Kinase case

% Column 1 (ancestral proteins)
% Column 2 (modern proteins)
Data(1).Name='4CSV_vs_Modern';
Data(1).TreePairTable={'4CSV'    'ABL1-H.sapien'
                       '4CSV'    'ABL2-H.sapien'
                       '4CSV'    'HCK-H.sapien'
                       '4CSV'    'LYN-H.sapien'
                       '4CSV'    'LCK-H.sapien'
                       '4CSV'    'SRC-H.sapien'
                       '4CSV'    'BTK-H.sapien'
                       '4CSV'    'BMX-H.sapien'
                       '4CSV'    'ITK-H.sapien'
                       '4CSV'    'HER1'};
Data(2).Name='4UEU_vs_Modern';    
Data(2).TreePairTable={'4UEU'    'ABL1-H.sapien'
                       '4UEU'    'ABL2-H.sapien'
                       '4UEU'    'HCK-H.sapien'
                       '4UEU'    'LYN-H.sapien'
                       '4UEU'    'LCK-H.sapien'
                       '4UEU'    'SRC-H.sapien'
                       '4UEU'    'BTK-H.sapien'
                       '4UEU'    'BMX-H.sapien'
                       '4UEU'    'ITK-H.sapien'
                       '4UEU'    'HER1'};

% Branch 73 is the node correspond to 4CSV and 4UEU in the published
% phylogenetic tree
Data(3).Name='Branch73_vs_Modern';
Data(3).TreePairTable={'Branch73'    'ABL1-H.sapien'
                       'Branch73'    'ABL2-H.sapien'
                       'Branch73'    'HCK-H.sapien'
                       'Branch73'    'LYN-H.sapien'
                       'Branch73'    'LCK-H.sapien'
                       'Branch73'    'SRC-H.sapien'
                       'Branch73'    'BTK-H.sapien'
                       'Branch73'    'BMX-H.sapien'
                       'Branch73'    'ITK-H.sapien'
                       'Branch73'    'HER1'};

Data(1).MSSATable={
'ABL1-H.sapien'  '_---__---EEEEEEEE______EEEEEETT----EEEEEE___S---S-----_S-THHHHHH-HHHHHHH___TTB__EEEEE_-SSS_-EEEEE__TTEEHHHHHH--_-TTTS--_HHHHHHHHHHHHHHHHHHHHTT---EEETT_SGGGEEEEG--------G----GEEEE_GGGEEE-_------TT_EE__T--T__-__TTTS_HHHHH-------H_EE_HHHHHHHHHHHHHHHHT-S__SSTT__HHHHHHHH--H-TT______T---T--__HHHH-HHHHHHT_SSGGGS__HHH--------'
'ABL2-H.sapien'  '_---__---EEEEEESGGGTTSSEEEEEEGG----EEEEEEEE_S---S-----_S-_HHHHHH-HHHHHTT___TTB__EEEEE_-SSSE-EEEEE__TT_BHHHHHH--_-TTTS--_HHHHHHHHHHHHHHHHHHHHTT---___S__SGGGEEE_G--------G----G_EEE_______-__-_---__SSB___--___-B_GGG__HHHHH-------H____HHHHHHHHHHHHHHHHT-S__SSTT__GGGHHHHH--H-HT______T---T--__HHHH-HHHHHHT_SSGGGS__HHH--------'
'HCK-H.sapien'   '-_--__---EEEEEEEE__SSEEEEEEEETT----EEEEEEEE_T---T-----SB-_HHHHHH-HHHHHTT___TTB__EEEEE_-SSS_-EEEE___TT_BHHHHHH----HHT_--_HHHHHHHHHHHHHHHHHHHHTT---___SS_SGGGEEE_T--------T----__EEE_STTGGG----_---_HHHHTT_--SSS-S_GGGS_HHHHH-------H____HHHHHHHHHHHHHHHHT-S__SSTT__HHHHHHHH--H-HT______T---T--S_HHHH-HHHHHHT_SSGGGS__HHH--------'
'LYN-H.sapien'   '_---__---B_______EETTEEEE__BBSS----_B_EEE___T---T-----SS-_SSSTTT-HHHHTTTT_BTTB___EEEE_-SSS_-EEEE___TTEEHHHHTT----HHHS--_HHHHHHHHHHHHHHHHHHHHTT---_B_S__SGGGEEE_T--------T----S_EEE___SS_B-__-_---________--___-__GGGS_HHHHH-------T___SHHHHHHHHHHHHHHHHT-S__SSTTS_TTTHHHHH--H-TT______S---S--S_HHHH-HHHHTTT_SSTTTS__HHH--------'
'LCK-H.sapien'   '_---__---EEEEEEEE__TTEEEEEEEETT----EEEEEEEE_T---T-----TS-_HHHHHH-HHHHHHT___TTB__EEEEE_-SSS_-EEEEE__TT_BHHHHTT----HHT_--_HHHHHHHHHHHHHHHHHHHHTT---___S__SGGGEEE_T--------T----S_EEE_______-__-_---________--___-___TTS_HHHHH-------H____HHHHHHHHHHHHHHHHT-T__SSTT__HHHHHHHH--H-HT______T---T--__HHHH-HHHHHHT_SSGGGS__HHH--------'
'SRC-H.sapien'   '-_--__---EEEEEEEEE_SSEEEEEEEETT----EEEEEEEE_T---T-----SS-_HHHHHH-HHHHHHH___TTB__EEEEE_-SSS_-EEEE___TTEEHHHHHS----HTT_--_HHHHHHHHHHHHHHHHHHHHTT---___S__SGGGEEE_G--------G----G_EEE__TTSTT--_-_---_HHHHTT_--STT-S_GGGS_HHHHH-------H____HHHHHHHHHHHHHHHTT-T__SSTT__HHHHHHHH--H-TT______T---T--__HHHH-HHHHHHT_SSGGGS__HHH--------'
'BTK-H.sapien'   '----__---EEEEEEEEE_SSEEEEEEEETT----EEEEEEEE_T---T-----_B-_HHHHHH-TTHHHHT___TTB__EEEEE_-SSS_-EEEEE__TT_BHHHHHH_---TS__--_HHHHHHHHHHHHHHHHHHHHTT---___S___GGGEEE_T--------T----__EEE_STTGGG----_---_HHHHSTT--STT-S_GGG__HHHHH-------H___SHHHHHHHHHHHHHHHHT-T__TTTTS_HHHHHHHH--T-TT______T---T--__HHHH-HHHHGGG_SSGGGS__HHH--------'
'BMX-H.sapien'   '_---__---EEEEEEEEEETTEEEEEEEETT----EEEEEEEE_B---T-----TB-_HHHHHH-HHHHHHH___TTB__EEEEE_-SSSE-EEEEE__TT_BHHHHHH----GGG_--_HHHHHHHHHHHHHHHHHHHHTT---EEESS_SGGGEEE_T--------T----__EEE__TT_EE--_-_---TT_EEE__--S__-__GGG__HHHHH-------HSEEETTHHHHHHHHHHHHHHT-T__TTTTS_HHHHHHHH--H-TT______T---T--S_HHHH-HHHHHTT_SSGGGS__HHH--------'
'ITK-H.sapien'   '_---__---EEEEEEEEEETTEEEEEEEETT----EEEEEEE__T---T-----SB-_HHHHHH-HHHHHHT___TTB__EEEEE_-SSS_-EEEEE__TT_BHHHHHH----TT__--_HHHHHHHHHHHHHHHHHHHHTT---___S___GGGEEE_G--------G----G_EEE__TTGGG----_---_HHHHSTT--STT-S_GGG__HHHHH-------H____HHHHHHHHHHHHHHHHT-S__TTTT__HHHHHHHH--H-TT______T---T--S_HHHH-HHHHHHT_SSGGGS__HHH--------'
'HER1'           '---------EEEEEEEEEETTEEEEEEEETT----EEEEEEE___---________-_HHHHHH-HHHHHHH___TTB__EEEEE______-EEEEE__TTEEHHHHHT----SS__--_HHHHHHHHHHHHHHHHHHHHSS--____S__SGGGEEESS__---__-S----__EEE_______--------________--___-__GGGS_HHHHH-------H___SHHHHHHHHHHHHHHHHH-___TTTTS_HHHHHHHH--H-S_______T---T--__HHHH-HHHHHHT_SSGGGS__HHH--------'
'Branch73'       '---------EEEEEEEEEESSEEEEEEEETT----EEEEEEEEET---T-----SS-HHHHHHH-HHHHHHHHEETTBEEEEEEEE-SSSE-EEEEEEETTEBHHHHHH----TTTS--EHHHHHHHHHHHHHHHHHHHHTT---EEESSESGGGEEEEG--------G----GEEEEESTTSTT--------HHHHHHHH--SHH-SEGGGSEHHHHH-------HEEEEHHHHHHHHHHHHHHHHT-SEESSTTSEHHHHHHHH--H-TTEEEEEET---T--SEHHHH-HHHHHHTESSGGGSEEHHH--------'
'4CSV'           '___B__GGGEEEEEEEE__SS__EEEEEETT--TTEEEEEEE__T---T-----SS-_HHHHHH-HHHHHTT___TTB__EEEEE_SSSS_-EEEEE__TT_BHHHHHTT__-____--_HHHHHHHHHHHHHHHHHHHHTT---___S___GGGEEE_G--------G----G_EEE_______-__-_---________--___-__GGGS_HHHHH-------H____HHHHHHHHHHHHHHHHTTS__SSTT__HHHHHHHH--H-TT______T---T--__HHHH-HHHHHHT_SSGGGS__HHHHHHHHHT_'
'4UEU'           '_TTB__GGGB___B__S__SS____EEEETT--TTEEEEEEE__S---S-----SS-_HHHHHH-HHHHHTT___TTB__EEEEE_SSSS_-EEEEE__TTEEHHHHHHSTT-GGG_--_HHHHHHHHHHHHHHHHHHHHTT---_B_S__SGGGEEE_G--------G----G_EEE___TT_B-__-_---________--___-__GGGS_HHHHH-------H____HHHHHHHHHHHHHHHHTTT__SSTT__SHHHHHHH--H-TT______T---T--S_HHHH-HHHHHHT_SSGGGS__HHHHHHHHHS_'};
Data(2).MSSATable=Data(1).MSSATable;
Data(3).MSSATable=Data(1).MSSATable;


% Enzyme case

% Column 1 (ancestral proteins)
% Column 2 (modern proteins)
Data(4).Name='4ZV1_vs_Modern';
Data(4).TreePairTable={'4ZV1' 'E1Enteroco'
                       '4ZV1' 'R1Streptoc'
                       '4ZV1' 'H1Escheric'
                       '4ZV1' 'Q1Escheric'
                       '4ZV1' 'Q1Enteroco'
                       '4ZV1' 'R1Lactococ'};

Data(5).Name='4ZV2_vs_Modern';
Data(5).TreePairTable={'4ZV2' 'E1Enteroco'
                       '4ZV2' 'R1Streptoc'
                       '4ZV2' 'H1Escheric'
                       '4ZV2' 'Q1Escheric'
                       '4ZV2' 'Q1Enteroco'
                       '4ZV2' 'R1Lactococ'};
% Branch 326 is the node correspond to 4ZV1 and 4ZV2 in the published
% phylogenetic tree
Data(6).Name='Branch326_vs_Modern';
Data(6).TreePairTable={'Branch326' 'E1Enteroco'
                       'Branch326' 'R1Streptoc'
                       'Branch326' 'H1Escheric'
                       'Branch326' 'Q1Escheric'
                       'Branch326' 'Q1Enteroco'
                       'Branch326' 'R1Lactococ'};
           
Data(4).MSSATable={
    'E1Enteroco' '------HTEEEEEE_SEETTTEEE--TTTEEESHHHHHHHHHHHHH-_-T_EEEEEE__TTTHHHHHHHTS_SEE_SS_B__TTGGGTEEE____EEE_EEEEEETT_-S__SGGG__T_EEEEETT_HHHHHHHHHS----T---EEEESSHHHHHHHHHTTS_SEEEEEHHHHHHHHTT_T-----------EES_BSS_EEE__EEETT_HHHHHHHHHHHHHHHHHTHHHHHHHHH_--'
    'R1Streptoc' '------HTEEEEEE_S_BTTTBEE--TEEEEESHHHHHHHHHHHHH---T_EEEEEE__HHHHHHHHHTTS_SEE_SS_B__HHHHTTEEE_S__EE__EEEEEEGGG-T_SSGGGG___EEEEETTSHHHHHHHHH_----S---EEEES_HHHHHHHHHTTSSSEEEEEHHHHHHHHHH_T-------_---_________EE__EEESS_HHHHHHHHHHHHHHHHTTHHHHHHHHHH--'
    'H1Escheric' '------_SS_EEEE___BTTTBEE_-TT__EESHHHHHHHHHHHHT---T___EEEE__HHHHHHHHHTTS_SEE_SS_B__HHHHTTSEE_S__EE__EEEEEESS_-__SSGGGGTT_EEEEETTSHHHHHHHHHT----G-_-EEEESSHHHHHHHHHHTSSSEEEEEHHHHHHHTTTSG----------____TTTS_SEE__EE_SS_HHHHHHHHHHHHHHHHTTHHHHHHHTT_--'
    'Q1Escheric' '------___EEEEEESSBTTTBEE--ETTEEESHHHHHHHHHHHHH---T__EEEEEE_GGGHHHHHHTTSSSEEEEEEE__HHHHTTSEE_S__EEEEEEEEEETT__S_SSSTTTTT_EEEEETTSHHHHHHHHH___--S---EEEESSHHHHHHHHHTTS_SEEEEEHHHHHHHHHTTT-----------EEEEEEEEEEEEEEE_TT_HHHHHHHHHHHHHHHHTSHHHHHHHHHH--'
    'Q1Enteroco' '------SS_EEEEE_S_BTTTBEE_-TTS_EESHHHHHHHHHHHHH---T__EEEEE__HHHHHHHHHHTSSSEE_SS_B__TTGGGTEEE_S__EEEEEEEEEETT__S_SSGGGGTT_EEEEETTSHHHHHHHHHH----H__-EEEESSHHHHHHHHHHTS_SEEEEEHHHHHHHHHTT_-----------__S__EEEEEE__EEETT_HHHHHHHHHHHHHHHHHTHHHHHHHHH___'
    'R1Lactococ' '------SS_EEEEE_SSBTTTBEE_-TTS_EE_HHHHHHHHHHHHH---T__EEEEE__HHHHHHHHHTT_SSEE_SS_B__TTGGGTEEE_S__EEEEEEEEEETT__S_SSGGGGTT_EEEEETTSHHHHHHHHHH----H-_-EEEES_HHHHHHHHHHTSSSEEEEEHHHHHHHHHHT_--------_--__S__EEEEEE__EEETT_HHHHHHHHHHHHHHHHHSHHHHHHHHHS__'
    'Branch326'  '------HHEEEEEEESEETTTEEE--SSTEEESHHHHHHHHHHHHH---TEEEEEEEEEHHHHHHHHHHTSHSEEESSEEEEHHHHHTEEEESEEEEEEEEEEEETTH-SESSTHHHHHEEEEEETTSHHHHHHHHHH----T---EEEESSHHHHHHHHHHTSHSEEEEEHHHHHHHHHHHH-----------EEEEEEEEEEEEEEEETTEHHHHHHHHHHHHHHHHHSHHHHHHHHHH--'
    '4ZV1'       '________EEEEEE_SSBTTTBEE_-TTS_EESHHHHHHHHHHHHH---T_EEEEEE__GGGHHHHHHTTS_SEE_SS_B__HHHHTTEEE_S__EEE_EEEEEETT__S_SSGGGGTT_EEEEETTSHHHHHHHHHHHH--HT_EEEEESSHHHHHHHHHHTS_SEEEEEHHHHHHHHHT_T------TS_EEEES____SEEE__EEETT_HHHHHHHHHHHHHHHHSSHHHHHHHHHH__'
    '4ZV2'       '_______SEEEEEE_SSBTTTBEE_-TTS_EESHHHHHHHHHHHHH---T_EEEEEE__GGGHHHHHHTTS_SEE_SS_B__TTGGGTSEE_S__EEE_EEEEEETT__S_SSGGGGTT_EEEEETTSHHHHHHHHH___--___EEEEESSHHHHHHHHHHTS_SEEEEEHHHHHHHHHH_G------GG_EEEES____SEEE__EE_TT_HHHHHHHHHHHHHHHHSSHHHHHHHHHH__'};
Data(5).MSSATable=Data(4).MSSATable;
Data(6).MSSATable=Data(4).MSSATable;



% DNA-binding domain case

% Column 1 (ancestral proteins)
% Column 2 (modern proteins)
Data(7).Name='5CBX_vs_Modern';
Data(7).TreePairTable={'5CBX' 'RatNovAR'
                       '5CBX' 'HomSapPR'
                       '5CBX' 'RatNovGR'
                       '5CBX' 'HomSapGR'
                       '5CBX' 'HomSapMR'
                       '5CBX' 'HomSapERa'
                       '5CBX' 'HomSapRXRa'};

Data(8).Name='5CBY_vs_Modern';
Data(8).TreePairTable={'5CBY' 'RatNovAR'
                       '5CBY' 'HomSapPR'
                       '5CBY' 'RatNovGR'
                       '5CBY' 'HomSapGR'
                       '5CBY' 'HomSapMR'
                       '5CBY' 'HomSapERa'
                       '5CBY' 'HomSapRXRa'};

Data(9).Name='5CBZ_vs_Modern';
Data(9).TreePairTable={'5CBZ' 'RatNovAR'
                       '5CBZ' 'HomSapPR'
                       '5CBZ' 'RatNovGR'
                       '5CBZ' 'HomSapGR'
                       '5CBZ' 'HomSapMR'
                       '5CBZ' 'HomSapERa'
                       '5CBZ' 'HomSapRXRa'};

Data(10).Name='Branch82_vs_Modern';
% Branch 82 is the node correspond to 5CBY in the published
% phylogenetic tree
Data(10).TreePairTable={'Branch82' 'RatNovAR'
                        'Branch82' 'HomSapPR'
                        'Branch82' 'RatNovGR'
                        'Branch82' 'HomSapGR'
                        'Branch82' 'HomSapMR'
                        'Branch82' 'HomSapERa'
                        'Branch82' 'HomSapRXRa'};

% Branch 83 is the node correspond to 5CBX in the published
% phylogenetic tree
Data(11).Name='Branch83_vs_Modern';
Data(11).TreePairTable={'Branch83' 'RatNovAR'
                        'Branch83' 'HomSapPR'
                        'Branch83' 'RatNovGR'
                        'Branch83' 'HomSapGR'
                        'Branch83' 'HomSapMR'
                        'Branch83' 'HomSapERa'
                        'Branch83' 'HomSapRXRa'};

% Branch 103 is the node correspond to 5CBZ in the published
% phylogenetic tree
Data(12).Name='Branch103_vs_Modern';    
Data(12).TreePairTable={'Branch103' 'RatNovAR'
                        'Branch103' 'HomSapPR'
                        'Branch103' 'RatNovGR'
                        'Branch103' 'HomSapGR'
                        'Branch103' 'HomSapMR'
                        'Branch103' 'HomSapERa'
                        'Branch103' 'HomSapRXRa'};
           

Data(7).MSSATable={ 
'RatNovAR'           '_____B_SSS_SBEEEEETTEEEEHHHHHHHHHHHTS______SSSS_____TTTTTT_HHHHHHHHHHHT__S______________SSSSS____B_SBSB__HHHHHHHHHHHHTT_____SSSS____STTTGGG_HHHHHHHHHHTSB_________'
'HomSapPR'           '--___B_TTT_SB__EEETTEEE_HHHHHHHHHHHHSS_____SS_S_____TTTTTT_HHHHHHHHHHHT_________---___B_TTTSSB__EEETTEEE_HHHHHHHHHHHHS______SS_S_____TTTTTT_HHHHHHHHHHHT____-____-'
'RatNovGR'           '_--__B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS_____TTGGGT_HHHHHHHHHHHT__TT-----__--__B_TTT_SB__EEETTEEE_HHHHHHHHHHHHHT_____SSSS____STTGGGT_HHHHHHHHHHTT__TT------'
'HomSapGR'           '_____B_TTT_SB__EEETTEEE_HHHHHHHHHHHTS______SSSS_____TTGGGT_HHHHHHHHHHHT_______________B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS____STTTTTT_HHHHHHHHHHHT__________'
'HomSapMR'           '_____B_TTT_SB__SEETTEE__HHHHHHHHHHHHS______SSSS_____TTTTTT_HHHHHHHHHHTT____------_____B_TTT_SB__SEETTEE__HHHHHHHHHHHHT______SSSS_____GGGTTT_HHHHHHHHHHHT____------'
'HomSapERa'          '_____B_TTT_SB__EEETTEEE_HHHHHHHHHHHSS______SSSS____STTTTTT_HHHHHHHHHHHT_______________B_TTT_SB__EEETTEEE_HHHHHHHHHHHT_______SSSS_____TTGGGT_HHHHHHHHHHTT__________'
'HomSapRXRa'         '__-SEE_TTT__EESEEETTEEE_HHHHHHHHHHHHTT_____SS_S____STTTTTS_HHHHHHHHHHTT__GG--_________B_TTT_SB__SEETTEE__HHHHHHHHHHHHTT_____________STTTGGG_HHHHHHHHHHTT__GG--____'
'Branch82'           '---EEBETTTESBEEEEETTEEEEHHHHHHHHHHHHSHEEEEESSSSEEEEHTTTTTTEHHHHHHHHHHHTEETH---------EEBETTTSSBEEEEETTEEEEHHHHHHHHHHHHSHEEEEESSSSEEEESTTTTTTEHHHHHHHHHHHTEEHH------'
'Branch83'           '---EEBETTTESBEEEEETTEEEEHHHHHHHHHHHHSHEEEEESSSSEEEEHTTTTTTEHHHHHHHHHHHTEETH---------EEBETTTSSBEEEEETTEEEEHHHHHHHHHHHHSHEEEEESSSSEEEESTTTTTTEHHHHHHHHHHHTEEHH------'
'Branch103'          '---EEBETTTESBEEEEETTEEEEHHHHHHHHHHHHSHEEEEESSSSEEEEHTTTTTTEHHHHHHHHHHHTEETH---------EEBETTTSSBEEEEETTEEEEHHHHHHHHHHHHSHEEEEESSSSEEEESTTTTTTEHHHHHHHHHHHTEEHH------'
'5CBX'             '_____B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS_____TTTTTT_HHHHHHHHHHTT__TT___________B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS____STTTTTT_HHHHHHHHHHTT__S_______'
'5CBY'             '_____B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS____STTTTTT_HHHHHHHHHHHT__TT___________B_TTT_SB__EEETTEEE_HHHHHHHHHHHHT______SSSS_____TTGGGT_HHHHHHHHHHTT__________'
'5CBZ'             '_____B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS____STTTTTT_HHHHHHHHHHTT__S____________B_TTT_SB__EEETTEEE_HHHHHHHHHHHHS______SSSS____STTTTTT_HHHHHHHHHHTT___S______'};

Data(8).MSSATable =Data(7).MSSATable;
Data(9).MSSATable =Data(7).MSSATable;
Data(10).MSSATable=Data(7).MSSATable;
Data(11).MSSATable=Data(7).MSSATable;
Data(12).MSSATable=Data(7).MSSATable;













