//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#ifndef CFEM_QUADRATURE_TPP
#define CFEM_QUADRATURE_TPP
#include "Quadrature.h"
#endif

namespace cfem
{

    inline QuadratureRule QuadratureFactory::getGaussLengendre1DRule(GaussLegendre1DRule rule) {
        switch(rule) {
            case GaussLegendre1DRule::GAUSS_1Pt:
                return create1DRuleGaussLegendre1();
            case GaussLegendre1DRule::GAUSS_2Pt:
                return create1DRuleGaussLegendre2();
            case GaussLegendre1DRule::GAUSS_3Pt:
                return create1DRuleGaussLegendre3();
            case GaussLegendre1DRule::GAUSS_4Pt:
                return create1DRuleGaussLegendre4();
            case GaussLegendre1DRule::GAUSS_5Pt:
                return create1DRuleGaussLegendre5();
            default:
                CFEM_ERROR << "Requested 1D quadrature not implemented. Supported: 1, 2, 3, 4, 5.";
        }
        return {};
    }

    inline QuadratureRule QuadratureFactory::getTriangleRule(TriangleRule rule)
    {
        switch (rule)
        {
        case TriangleRule::GAUSS_1Pt:
            return createTriangleRuleGauss1();
        case TriangleRule::GAUSS_3Pt:
            return createTriangleRuleGauss3();
        case TriangleRule::GAUSS_4Pt:
            return createTriangleRuleGauss4();
        case TriangleRule::GAUSS_6Pt:
            return createTriangleRuleGauss6();
        case TriangleRule::GAUSS_7Pt:
            return createTriangleRuleGauss7();
        case TriangleRule::GAUSS_12Pt:
            return createTriangleRuleGauss12();
        case TriangleRule::GAUSS_16Pt:
            return createTriangleRuleGauss16();
        case TriangleRule::GAUSS_19Pt:
            return createTriangleRuleGauss19();
        case TriangleRule::GAUSS_25Pt:
            return createTriangleRuleGauss25();
        default:
            CFEM_ERROR << "Requested triangle quadrature not implemented. Supported: 1, 3, 4, 6, 7, 12, 16, 19.";
        }
        return {};
    }

    inline QuadratureRule QuadratureFactory::getQuadRule(CfemInt nPoints1D) {
        QuadratureRule rule;
        QuadratureRule base;

        switch (nPoints1D)
        {
        case 1:
            base = create1DRuleGaussLegendre1(); break;
        case 4:
            base = create1DRuleGaussLegendre2(); break;
        case 9:
            base = create1DRuleGaussLegendre3(); break;
        case 16:
            base = create1DRuleGaussLegendre4(); break;
        case 25:
            base = create1DRuleGaussLegendre5(); break;
        
        default:
            CFEM_ERROR << "Requested Quad quadrature not implemented. Supported: 1, 4, 9, 16, 25";
        }
        // Create the Hex rule
        for (const auto& px : base.getPoints()) {
            for (const auto& py : base.getPoints()) {
                rule.m_points.push_back({
                    {px.xi.x, py.xi.x, 0.0},
                    px.weight * py.weight 
                });
            }
        }
        return rule;
    }

    inline QuadratureRule QuadratureFactory::getTetRule(TetRule rule)
    {
        switch (rule)
        {
        case TetRule::GAUSS_1Pt:
            return createTetRuleGauss1();
        case TetRule::GAUSS_4Pt:
            return createTetRuleGauss4_deg2();
        case TetRule::GAUSS_5Pt:
            return createTetRuleGauss5_deg3();
        case TetRule::GAUSS_11Pt:
            return createTetRuleGauss11_deg4();
        case TetRule::GAUSS_15Pt:
            return createTetRuleGauss15_deg5();
        case TetRule::GAUSS_29Pt:
            return createTetRuleGauss29_deg6();
        case TetRule::GAUSS_33Pt:
            return createTetRuleGauss33_deg7();
        default:
            CFEM_ERROR << "Requested tetra quadrature not implemented. Supported: 1, 4, 5, 11, 15, 29, 33.";
        }
        return {};
    }

    /** @brief Produit tensoriel entre un Triangle (2D) et une ligne (1D) */
    inline QuadratureRule QuadratureFactory::getPrismRule(TriangleRule triRule, GaussLegendre1DRule element1DRule) {
        QuadratureRule rule;
        QuadratureRule baseTri = getTriangleRule(triRule);
        QuadratureRule baseLine = getGaussLengendre1DRule(element1DRule);

        for (const auto& pTri : baseTri.getPoints()) {
            for (const auto& pLine : baseLine.getPoints()) {
                rule.m_points.push_back({
                    {pTri.xi.x, pTri.xi.y, pLine.xi.x}, // xi, eta du triangle, zeta de la ligne
                    pTri.weight * pLine.weight
                });
            }
        }
        return rule;
    }

    inline QuadratureRule QuadratureFactory::getHexRule(CfemInt nPoints1D) {
        QuadratureRule rule;
        QuadratureRule base;

        switch (nPoints1D)
        {
        case 1:
            base = create1DRuleGaussLegendre1(); break;
        case 8:
            base = create1DRuleGaussLegendre2(); break;
        case 27:
            base = create1DRuleGaussLegendre3(); break;
        case 64:
            base = create1DRuleGaussLegendre4(); break;
        case 125:
            base = create1DRuleGaussLegendre5(); break;
        
        default:
            CFEM_ERROR << "Requested Hex quadrature not implemented. Supported: 1, 8, 27, 64, 125";
        }

        for (const auto& px : base.getPoints()) {
            for (const auto& py : base.getPoints()) {
                for (const auto& pz : base.getPoints()) {
                    rule.m_points.push_back({
                        {px.xi.x, py.xi.x, pz.xi.x}, 
                        px.weight * py.weight * pz.weight
                    });
                }
            }
        }
        return rule;
    }


    inline QuadratureRule QuadratureFactory::create1DRuleGaussLegendre1(){
        QuadratureRule rule;
        rule.m_points.push_back({{0, 0, 0}, 2.0});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::create1DRuleGaussLegendre2(){
        QuadratureRule rule;
        const CfemReal x = 0.57735026918962576451L; // 1/sqrt(3)
        rule.m_points.push_back({{-x, 0, 0}, 1.0});
        rule.m_points.push_back({{ x, 0, 0}, 1.0});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::create1DRuleGaussLegendre3(){
        QuadratureRule rule;
        const CfemReal x = 0.77459666924148337704L; // sqrt(0.6)
        rule.m_points.push_back({{-x, 0, 0}, 5.0L/9.0L});
        rule.m_points.push_back({{ 0, 0, 0}, 8.0L/9.0});
        rule.m_points.push_back({{ x, 0, 0}, 5.0L/9.0L});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::create1DRuleGaussLegendre4(){
        QuadratureRule rule;
        const CfemReal x1 = 0.86113631159405257522L;
        const CfemReal x2 = 0.33998104358485626480L;
        const CfemReal w1 = 0.34785484513745385737L;
        const CfemReal w2 = 0.65214515486254614263L;
        rule.m_points.push_back({{-x1, 0, 0}, w1});
        rule.m_points.push_back({{-x2, 0, 0}, w2});
        rule.m_points.push_back({{ x2, 0, 0}, w2});
        rule.m_points.push_back({{ x1, 0, 0}, w1});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::create1DRuleGaussLegendre5(){
        QuadratureRule rule;
        const CfemReal x1 = 0.90617984593866399280L;  // 1.0/21.0*std::sqrt(245.0+14.0*std::sqrt(70.0));
        const CfemReal x2 = 0.53846931010568309104L;  // 1.0/21.0*std::sqrt(245.0-14.0*std::sqrt(70.0));
        const CfemReal w1 = 0.23692688505618908751L;  // 1.0/900.0*(322.0-13.0*std::sqrt(70.0));
        const CfemReal w2 = 0.47862867049936646804L;  // 1.0/900.0*(322.0+13.0*std::sqrt(70.0));
        const CfemReal w3 = 0.56888888888888888889L;  // 128.0/225.0;
        rule.m_points.push_back({{-x1, 0, 0}, w1});
        rule.m_points.push_back({{-x2, 0, 0}, w2});
        rule.m_points.push_back({{ 0,  0, 0}, w3});
        rule.m_points.push_back({{ x2, 0, 0}, w2});
        rule.m_points.push_back({{ x1, 0, 0}, w1});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss1()
    {
        QuadratureRule rule;
        const CfemReal a = 1.0L / 3.0L;
        rule.m_points.push_back({{a, a, 0.0}, 0.5});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss3()
    {
        QuadratureRule rule;
        constexpr CfemReal w = 1.0 / 6.0;
        rule.m_points.push_back({{1.0 / 6.0, 1.0 / 6.0, 0.0}, w});
        rule.m_points.push_back({{2.0 / 3.0, 1.0 / 6.0, 0.0}, w});
        rule.m_points.push_back({{1.0 / 6.0, 2.0 / 3.0, 0.0}, w});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss4()
    {
        QuadratureRule rule;
        constexpr CfemReal a = 1.0 / 3.0;
        constexpr CfemReal b = 1.0 / 5.0;
        constexpr CfemReal c = 3.0 / 5.0;
        constexpr CfemReal w1 = -27.0 / 96.0;
        constexpr CfemReal w2 = 25.0 / 96.0;

        rule.m_points.push_back({{a, a, 0.0}, w1});
        rule.m_points.push_back({{b, c, 0.0}, w2});
        rule.m_points.push_back({{b, b, 0.0}, w2});
        rule.m_points.push_back({{c, b, 0.0}, w2});
        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss6()
    {
        QuadratureRule rule;
        const CfemReal x1 = 0.10810301816807000669L;
        const CfemReal y1 = 0.44594849091596499665L;
        const CfemReal w1 = 0.11169079483900552733L;

        // Orbite 2
        const CfemReal x2 = 0.81684757298045851052L;
        const CfemReal y2 = 0.09157621350977074474L;
        const CfemReal w2 = 0.05497587182766113934L;

        rule.m_points = {
            {{x1, y1, 0.0}, w1},
            {{y1, x1, 0.0}, w1},
            {{y1, y1, 0.0}, w1}, // Ici xi + eta + lambda3 = 1 -> y1 + y1 + x1 = 1
            {{x2, y2, 0.0}, w2},
            {{y2, x2, 0.0}, w2},
            {{y2, y2, 0.0}, w2}
        };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss7()
    {
        QuadratureRule rule;
        const CfemReal a6 = 0.470142064105115089770441209513447L;
        const CfemReal b6 = 0.101286507323456338800987361915123L;
        const CfemReal c6 = 0.0661970763942530903688246939165759L;
        const CfemReal d6 = 0.06296959027241357629784197275009064183802671L;

        rule.m_points = {
            {{b6,         b6, 0.0}, d6},
            {{1.L-2.L*b6, b6, 0.0}, d6},
            {{b6,   1.-2.*b6, 0.0}, d6},
            {{a6,   1.-2.*a6, 0.0}, c6},
            {{a6,         a6, 0.0}, c6},
            {{1.-2.*a6,   a6, 0.0}, c6},
            {{1./3.,   1./3., 0.0}, 9./80.}
        };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss12()
    {
        QuadratureRule rule;
        rule.m_points = {
            {{6.30890144915022283403316028708191e-02,6.30890144915022283403316028708191e-02, 0.0}, 2.54224531851034084604684045534344e-02},
            {{8.73821971016995543319336794258362e-01,6.30890144915022283403316028708191e-02, 0.0}, 2.54224531851034084604684045534344e-02},
            {{6.30890144915022283403316028708191e-02,8.73821971016995543319336794258362e-01, 0.0}, 2.54224531851034084604684045534344e-02},
            {{2.49286745170910421291638553107019e-01,2.49286745170910421291638553107019e-01, 0.0}, 5.83931378631896830126448056927897e-02},
            {{5.01426509658179157416722893785962e-01,2.49286745170910421291638553107019e-01, 0.0}, 5.83931378631896830126448056927897e-02},
            {{2.49286745170910421291638553107019e-01,5.01426509658179157416722893785962e-01, 0.0}, 5.83931378631896830126448056927897e-02},
            {{3.10352451033784405416607733956552e-01,5.31450498448169473532496716313981e-02, 0.0}, 4.14255378091867875967767282102212e-02},
            {{5.31450498448169473532496716313981e-02,3.10352451033784405416607733956552e-01, 0.0}, 4.14255378091867875967767282102212e-02},
            {{6.36502499121398647230142594412050e-01,3.10352451033784405416607733956552e-01, 0.0}, 4.14255378091867875967767282102212e-02},
            {{6.36502499121398647230142594412050e-01,5.31450498448169473532496716313981e-02, 0.0}, 4.14255378091867875967767282102212e-02},
            {{3.10352451033784405416607733956552e-01,6.36502499121398647230142594412050e-01, 0.0}, 4.14255378091867875967767282102212e-02},
            {{5.31450498448169473532496716313981e-02,6.36502499121398647230142594412050e-01, 0.0}, 4.14255378091867875967767282102212e-02}
        };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss16()
    {
        QuadratureRule rule;
        
        const CfemReal a10 = 0.081414823414553687942368971011662;
        const CfemReal b10 = 0.459292588292723156028815514494169;
        const CfemReal c10 = 0.658861384496479586755412997017072;
        const CfemReal d10 = 0.170569307751760206622293501491464;
        const CfemReal e10 = 0.898905543365938049083152898806802;
        const CfemReal f10 = 0.0505472283170309754584235505965989;
        const CfemReal g10 = 0.008394777409957605337213834539296;
        const CfemReal h10 = 0.263112829634638113421785786284643;
        const CfemReal i10 = 0.728492392955404281241000379176061;
        const CfemReal j10 = 0.0721578038388935841255455552445323;
        const CfemReal k10 = 0.0475458171336423123969480521942921;
        const CfemReal l10 = 0.0516086852673591251408957751460645;
        const CfemReal m10 = 0.0162292488115990401554629641708902;
        const CfemReal n10 = 0.01361515708721749713242234503695444183802671;

        rule.m_points =  {
            { {1.L/3.L, 1.L/3.L, 0.0}, j10},
    	    { {a10,  b10,   0.0}, k10},
            { {b10,  a10,   0.0}, k10},
            { {b10,  b10,   0.0}, k10},
            { {c10,  d10,   0.0}, l10},
            { {d10,  c10,   0.0}, l10},
            { {d10,  d10,   0.0}, l10},
            { {e10,  f10,   0.0}, m10},
            { {f10,  e10,   0.0}, m10},
            { {f10,  f10,   0.0}, m10},
            { {g10,  h10,   0.0}, n10},
            { {g10,  i10,   0.0}, n10},
            { {h10,  g10,   0.0}, n10},
            { {h10,  i10,   0.0}, n10},
            { {i10,  g10,   0.0}, n10},
            { {i10,  h10,   0.0}, n10}
        };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss19()
    {
        QuadratureRule rule;
        rule.m_points = {
            {{3.33333333333333333333333333333333e-01,3.33333333333333333333333333333333e-01, 0.0}, 4.85678981413994169096209912536443e-02},
            {{4.89682519198737627783706924836192e-01,4.89682519198737627783706924836192e-01, 0.0}, 1.56673501135695352684274156436046e-02},
            {{2.06349616025247444325861503276160e-02,4.89682519198737627783706924836192e-01, 0.0}, 1.56673501135695352684274156436046e-02},
            {{4.89682519198737627783706924836192e-01,2.06349616025247444325861503276160e-02, 0.0}, 1.56673501135695352684274156436046e-02},
            {{4.37089591492936637269930364435354e-01,4.37089591492936637269930364435354e-01, 0.0}, 3.89137705023871396583696781497019e-02},
            {{1.25820817014126725460139271129292e-01,4.37089591492936637269930364435354e-01, 0.0}, 3.89137705023871396583696781497019e-02},
            {{4.37089591492936637269930364435354e-01,1.25820817014126725460139271129292e-01, 0.0}, 3.89137705023871396583696781497019e-02},
            {{1.88203535619032730240961280467335e-01,1.88203535619032730240961280467335e-01, 0.0}, 3.98238694636051265164458871320226e-02},
            {{6.23592928761934539518077439065330e-01,1.88203535619032730240961280467335e-01, 0.0}, 3.98238694636051265164458871320226e-02},
            {{1.88203535619032730240961280467335e-01,6.23592928761934539518077439065330e-01, 0.0}, 3.98238694636051265164458871320226e-02},
    	    {{4.47295133944527098651065899662763e-02,4.47295133944527098651065899662763e-02, 0.0}, 1.27888378293490156308393992794999e-02},
    	    {{9.10540973211094580269786820067447e-01,4.47295133944527098651065899662763e-02, 0.0}, 1.27888378293490156308393992794999e-02},
            {{4.47295133944527098651065899662763e-02,9.10540973211094580269786820067447e-01, 0.0}, 1.27888378293490156308393992794999e-02},
            {{7.41198598784498020690079873523423e-01,3.68384120547362836348175987833851e-02, 0.0}, 2.16417696886446886446886446886446e-02},
    	    {{3.68384120547362836348175987833851e-02,7.41198598784498020690079873523423e-01, 0.0}, 2.16417696886446886446886446886446e-02},
    	    {{2.21962989160765695675102527693192e-01,7.41198598784498020690079873523423e-01, 0.0}, 2.16417696886446886446886446886446e-02},
            {{2.21962989160765695675102527693192e-01,3.68384120547362836348175987833851e-02, 0.0}, 2.16417696886446886446886446886446e-02},
            {{7.41198598784498020690079873523423e-01,2.21962989160765695675102527693192e-01, 0.0}, 2.16417696886446886446886446886446e-02},
            {{3.68384120547362836348175987833851e-02,2.21962989160765695675102527693192e-01, 0.0}, 2.16417696886446886446886446886446e-02}
        };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTriangleRuleGauss25()
    {
        QuadratureRule rule;
        rule.m_points = {
            {{3.33333333333333333333333333333333e-01,3.33333333333333333333333333333333e-01, 0.0},  3.99472523706198539156235226066932e-02},
            {{4.25086210602090572969529511638044e-01,4.25086210602090572969529511638044e-01, 0.0},  3.55619011161886673196456436993290e-02},
            {{1.49827578795818854060940976723912e-01,4.25086210602090572969529511638044e-01, 0.0},  3.55619011161886673196456436993290e-02},
            {{4.25086210602090572969529511638044e-01,1.49827578795818854060940976723912e-01, 0.0},  3.55619011161886673196456436993290e-02},
    	    {{2.33088675100001907144663868959796e-02,2.33088675100001907144663868959796e-02, 0.0},  4.11190934523209775932331018123594e-03},
            {{9.53382264979999618571067226208041e-01,2.33088675100001907144663868959796e-02, 0.0},  4.11190934523209775932331018123594e-03},
            {{2.33088675100001907144663868959796e-02,9.53382264979999618571067226208041e-01, 0.0},  4.11190934523209775932331018123594e-03},
            {{6.28307400213492556420837666078834e-01,2.23766973576973006225686490268204e-01, 0.0},  2.27152961480850090035368146219665e-02},
            {{2.23766973576973006225686490268204e-01,6.28307400213492556420837666078834e-01, 0.0},  2.27152961480850090035368146219665e-02},
    	    {{1.47925626209534437353475843652962e-01,6.28307400213492556420837666078834e-01, 0.0},  2.27152961480850090035368146219665e-02},
    	    {{1.47925626209534437353475843652962e-01,2.23766973576973006225686490268204e-01, 0.0},  2.27152961480850090035368146219665e-02},
            {{6.28307400213492556420837666078834e-01,1.47925626209534437353475843652962e-01, 0.0},  2.27152961480850090035368146219665e-02},
            {{2.23766973576973006225686490268204e-01,1.47925626209534437353475843652962e-01, 0.0},  2.27152961480850090035368146219665e-02},
            {{6.11313826181397648918755002253901e-01,3.58740141864431464578155300723852e-01, 0.0},  1.86799281171526384131182495009876e-02},
            {{3.58740141864431464578155300723852e-01,6.11313826181397648918755002253901e-01, 0.0},  1.86799281171526384131182495009876e-02},
            {{2.99460319541708865030896970222470e-02,6.11313826181397648918755002253901e-01, 0.0},  1.86799281171526384131182495009876e-02},
            {{2.99460319541708865030896970222470e-02,3.58740141864431464578155300723852e-01, 0.0},  1.86799281171526384131182495009876e-02},
            {{6.11313826181397648918755002253901e-01,2.99460319541708865030896970222470e-02, 0.0},  1.86799281171526384131182495009876e-02},
            {{3.58740141864431464578155300723852e-01,2.99460319541708865030896970222470e-02, 0.0},  1.86799281171526384131182495009876e-02},
            {{8.21072069985629373373544413472177e-01,1.43295370426867145305856630617323e-01, 0.0},  1.54433284422819943912565385023143e-02},
            {{1.43295370426867145305856630617323e-01,8.21072069985629373373544413472177e-01, 0.0},  1.54433284422819943912565385023143e-02},
            {{3.56325595875034813205989559105000e-02,8.21072069985629373373544413472177e-01, 0.0},  1.54433284422819943912565385023143e-02},
            {{3.56325595875034813205989559105000e-02,1.43295370426867145305856630617323e-01, 0.0},  1.54433284422819943912565385023143e-02},
            {{8.21072069985629373373544413472177e-01,3.56325595875034813205989559105000e-02, 0.0},  1.54433284422819943912565385023143e-02},
            {{1.43295370426867145305856630617323e-01,3.56325595875034813205989559105000e-02, 0.0},  1.54433284422819943912565385023143e-02}
        };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss1()
    {
        QuadratureRule rule;
        rule.m_points =  { {{0.25, 0.25, 0.25}, 1.0 / 6.0} };
        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss4_deg2()
    {
        QuadratureRule rule;
        // Parameters for degree 2 / 4-point rule
        const CfemReal a = 0.58541019662496845446L; // (5.0L + 3.0L * sqrt5) / 20.0L;
        const CfemReal b = 0.13819660112501051518L; // (5.0L - sqrt5) / 20.0L; or 0.138196601125010515179541316563436
        const CfemReal w = 1.0 / 24.0;              // Each point carries 1/4 of total volume
        
        rule.m_points.resize(4);
        rule.m_points[0] = {{a, b, b}, w};
        rule.m_points[1] = {{b, a, b}, w};
        rule.m_points[2] = {{b, b, a}, w};
        rule.m_points[3] = {{b, b, b}, w};
        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss5_deg3()
    {
        QuadratureRule rule;
        const CfemReal w_center = -2.0L / 15.0L;
        const CfemReal w_outer = 3.0L / 40.0L;
        const CfemReal a = 0.5L;
        const CfemReal b = 1.0L / 6.0L;

        rule.m_points.resize(5);
        rule.m_points[0] = {{0.25L, 0.25L, 0.25L}, w_center};
        rule.m_points[1] = {{a, b, b}, w_outer};
        rule.m_points[2] = {{b, a, b}, w_outer};
        rule.m_points[3] = {{b, b, a}, w_outer};
        rule.m_points[4] = {{b, b, b}, w_outer};

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss11_deg4()
    {
        QuadratureRule rule;
        const CfemReal lCentre_1     = 0.25;

        const CfemReal lSommets_L123 =  1./14.;
        const CfemReal lSommets_L4   = 11./14.;

        const CfemReal lAretes_L12   = .39940357616679920499610214746164062310585988737476; //0.25 * (1. + sqrt(5./14.));
        const CfemReal lAretes_L34   = .10059642383320079500389785253835937689414011262524; //0.25 * (1. - sqrt(5./14.));

        const CfemReal lPoidsCentre  = -74. /  5625.;
        const CfemReal lPoidsSommets = 343. / 45000.;
        const CfemReal lPoidsAretes  =  56. /  2250.;


        rule.m_points.resize(11);
        rule.m_points[0] = {{lCentre_1, lCentre_1, lCentre_1}, lPoidsCentre};

        rule.m_points[1] = {{lSommets_L123, lSommets_L123, lSommets_L123}, lPoidsSommets};
        rule.m_points[2] = {{lSommets_L4,   lSommets_L123, lSommets_L123}, lPoidsSommets};
        rule.m_points[3] = {{lSommets_L123, lSommets_L4,   lSommets_L123}, lPoidsSommets};
        rule.m_points[4] = {{lSommets_L123, lSommets_L123, lSommets_L4}, lPoidsSommets};

        rule.m_points[5] = {{lAretes_L12, lAretes_L12, lAretes_L34}, lPoidsAretes};
        rule.m_points[6] = {{lAretes_L12, lAretes_L34, lAretes_L12}, lPoidsAretes};
        rule.m_points[7] = {{lAretes_L34, lAretes_L12, lAretes_L12}, lPoidsAretes};

        rule.m_points[8] = {{lAretes_L34, lAretes_L12, lAretes_L34}, lPoidsAretes};
        rule.m_points[9] = {{lAretes_L12, lAretes_L34, lAretes_L34}, lPoidsAretes};
        rule.m_points[10]= {{lAretes_L34, lAretes_L34, lAretes_L12}, lPoidsAretes};

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss15_deg5()
    {
        QuadratureRule rule;

        const CfemReal a4 = .25;
        const CfemReal b4 = .91971078052723032788845135300e-1;
        const CfemReal bb4 = .3197936278296299083876254529344183802671;
        const CfemReal c4 = .724086765841830901633464594098;
        const CfemReal cc4 = .40619116511110274837123641196e-1;
        const CfemReal d4 = .56350832689629155741036730011e-1;
        const CfemReal e4 = .443649167310370844258963269989;

        const CfemReal f4 = .197530864197530864197530864198e-1;
        const CfemReal g4 = .119895139631697700000000000000e-1;
        const CfemReal gg4 = .115113678710453975468333333334e-1;
        const CfemReal h4 = .881834215167548500881834215168e-2;

        rule.m_points.resize(15);
        rule.m_points[0] =  {{a4, a4, a4}, f4};
        rule.m_points[1] =  {{b4, b4, b4}, g4};
        rule.m_points[2] =  {{b4, b4, c4}, g4};
        rule.m_points[3] =  {{b4, c4, b4}, g4};
        rule.m_points[4] =  {{c4, b4, b4}, g4};
        rule.m_points[5] =  {{bb4, bb4, bb4}, gg4};
        rule.m_points[6] =  {{bb4, bb4, cc4}, gg4};
        rule.m_points[7] =  {{bb4, cc4, bb4}, gg4};
        rule.m_points[8] =  {{cc4, bb4, bb4}, gg4};
        rule.m_points[9] =  {{d4, d4, e4}, h4};
        rule.m_points[10] = {{d4, e4, d4}, h4};
        rule.m_points[11] = {{e4, d4, d4}, h4};
        rule.m_points[12] = {{d4, e4, e4}, h4};
        rule.m_points[13] = {{e4, d4, e4}, h4};
        rule.m_points[14] = {{e4, e4, d4}, h4};

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss29_deg6()
    {
        QuadratureRule rule;
        rule.m_points.resize(29);

        rule.m_points[0] =  {{2.50000000000000000000000000000000e-01,2.50000000000000000000000000000000e-01,2.50000000000000000000000000000000e-01}, 1.50668817433579497383277309990912e-02 };
        rule.m_points[1] =  {{8.27719248047929541260063824577531e-01,5.74269173173568195799787251408230e-02,5.74269173173568195799787251408230e-02}, 3.18663904649853147632014415654494e-03 };
        rule.m_points[2] =  {{5.74269173173568195799787251408230e-02,8.27719248047929541260063824577531e-01,5.74269173173568195799787251408230e-02}, 3.18663904649853147632014415654494e-03 };
        rule.m_points[3] =  {{5.74269173173568195799787251408230e-02,5.74269173173568195799787251408230e-02,8.27719248047929541260063824577531e-01}, 3.18663904649853147632014415654494e-03 };
        rule.m_points[4] =  {{5.74269173173568195799787251408230e-02,5.74269173173568195799787251408230e-02,5.74269173173568195799787251408230e-02}, 3.18663904649853147632014415654494e-03 };
        rule.m_points[5] =  {{5.13518841255633944430944072477870e-02,4.86051028570607278709198710768507e-01,2.31298543651914663423853440991853e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[6] =  {{5.13518841255633944430944072477870e-02,2.31298543651914663423853440991853e-01,4.86051028570607278709198710768507e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[7] =  {{5.13518841255633944430944072477870e-02,2.31298543651914663423853440991853e-01,2.31298543651914663423853440991853e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[8] =  {{4.86051028570607278709198710768507e-01,5.13518841255633944430944072477870e-02,2.31298543651914663423853440991853e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[9] =  {{4.86051028570607278709198710768507e-01,2.31298543651914663423853440991853e-01,5.13518841255633944430944072477870e-02}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[10] = {{4.86051028570607278709198710768507e-01,2.31298543651914663423853440991853e-01,2.31298543651914663423853440991853e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[11] = {{2.31298543651914663423853440991853e-01,5.13518841255633944430944072477870e-02,4.86051028570607278709198710768507e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[12] = {{2.31298543651914663423853440991853e-01,5.13518841255633944430944072477870e-02,2.31298543651914663423853440991853e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[13] = {{2.31298543651914663423853440991853e-01,4.86051028570607278709198710768507e-01,5.13518841255633944430944072477870e-02}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[14] = {{2.31298543651914663423853440991853e-01,4.86051028570607278709198710768507e-01,2.31298543651914663423853440991853e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[15] = {{2.31298543651914663423853440991853e-01,2.31298543651914663423853440991853e-01,5.13518841255633944430944072477870e-02}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[16] = {{2.31298543651914663423853440991853e-01,2.31298543651914663423853440991853e-01,4.86051028570607278709198710768507e-01}, 7.26915640111093824271522019500777e-03 };
        rule.m_points[17] = {{2.96753812969025994633566939960374e-01,6.08107989401528086074390221633350e-01,4.75690988147229596460214192031380e-02}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[18] = {{2.96753812969025994633566939960374e-01,4.75690988147229596460214192031380e-02,6.08107989401528086074390221633350e-01}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[19] = {{2.96753812969025994633566939960374e-01,4.75690988147229596460214192031380e-02,4.75690988147229596460214192031380e-02}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[20] = {{6.08107989401528086074390221633350e-01,2.96753812969025994633566939960374e-01,4.75690988147229596460214192031380e-02}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[21] = {{6.08107989401528086074390221633350e-01,4.75690988147229596460214192031380e-02,2.96753812969025994633566939960374e-01}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[22] = {{6.08107989401528086074390221633350e-01,4.75690988147229596460214192031380e-02,4.75690988147229596460214192031380e-02}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[23] = {{4.75690988147229596460214192031380e-02,2.96753812969025994633566939960374e-01,6.08107989401528086074390221633350e-01}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[24] = {{4.75690988147229596460214192031380e-02,2.96753812969025994633566939960374e-01,4.75690988147229596460214192031380e-02}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[25] = {{4.75690988147229596460214192031380e-02,6.08107989401528086074390221633350e-01,2.96753812969025994633566939960374e-01}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[26] = {{4.75690988147229596460214192031380e-02,6.08107989401528086074390221633350e-01,4.75690988147229596460214192031380e-02}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[27] = {{4.75690988147229596460214192031380e-02,4.75690988147229596460214192031380e-02,2.96753812969025994633566939960374e-01}, 4.30194599366527767587297639177519e-03 };
        rule.m_points[28] = {{4.75690988147229596460214192031380e-02,4.75690988147229596460214192031380e-02,6.08107989401528086074390221633350e-01}, 4.30194599366527767587297639177519e-03 };

        return rule;
    }

    inline QuadratureRule QuadratureFactory::createTetRuleGauss33_deg7()
    {
        QuadratureRule rule;
        rule.m_points.resize(33);

        rule.m_points[0] =  {{2.50000000000000000000000000000000e-01,2.50000000000000000000000000000000e-01,2.50000000000000000000000000000000e-01},-1.04010315215862696536957218201166e-01};
        rule.m_points[1] =  {{3.44546088363740739314173508292220e-01,2.18484637212086420228608830569260e-01,2.18484637212086420228608830569260e-01}, 3.78340645050332670960907661584800e-02};
        rule.m_points[2] =  {{2.18484637212086420228608830569260e-01,3.44546088363740739314173508292220e-01,2.18484637212086420228608830569260e-01}, 3.78340645050332670960907661584800e-02};
        rule.m_points[3] =  {{2.18484637212086420228608830569260e-01,2.18484637212086420228608830569260e-01,3.44546088363740739314173508292220e-01}, 3.78340645050332670960907661584800e-02};
        rule.m_points[4] =  {{2.18484637212086420228608830569260e-01,2.18484637212086420228608830569260e-01,2.18484637212086420228608830569260e-01}, 3.78340645050332670960907661584800e-02};
        rule.m_points[5] =  {{8.19628940628295702134923655311801e-01,6.01236864572347659550254482293996e-02,6.01236864572347659550254482293996e-02}, 3.48619089127207855486883209460898e-03};
        rule.m_points[6] =  {{6.01236864572347659550254482293996e-02,8.19628940628295702134923655311801e-01,6.01236864572347659550254482293996e-02}, 3.48619089127207855486883209460898e-03};
        rule.m_points[7] =  {{6.01236864572347659550254482293996e-02,6.01236864572347659550254482293996e-02,8.19628940628295702134923655311801e-01}, 3.48619089127207855486883209460898e-03};
        rule.m_points[8] =  {{6.01236864572347659550254482293996e-02,6.01236864572347659550254482293996e-02,6.01236864572347659550254482293996e-02}, 3.48619089127207855486883209460898e-03};
        rule.m_points[9] =  {{2.96366112714042992044334920184189e-01,6.03031636054078023896939525318169e-01,5.03011256159394920293627772488211e-02}, 4.66737035605682874228409536944228e-03};
        rule.m_points[10] = {{2.96366112714042992044334920184189e-01,5.03011256159394920293627772488211e-02,6.03031636054078023896939525318169e-01}, 4.66737035605682874228409536944228e-03};
        rule.m_points[11] = {{2.96366112714042992044334920184189e-01,5.03011256159394920293627772488211e-02,5.03011256159394920293627772488211e-02}, 4.66737035605682874228409536944228e-03};
        rule.m_points[12] = {{6.03031636054078023896939525318169e-01,2.96366112714042992044334920184189e-01,5.03011256159394920293627772488211e-02}, 4.66737035605682874228409536944228e-03};
        rule.m_points[13] = {{6.03031636054078023896939525318169e-01,5.03011256159394920293627772488211e-02,2.96366112714042992044334920184189e-01}, 4.66737035605682874228409536944228e-03};
        rule.m_points[14] = {{6.03031636054078023896939525318169e-01,5.03011256159394920293627772488211e-02,5.03011256159394920293627772488211e-02}, 4.66737035605682874228409536944228e-03};
        rule.m_points[15] = {{5.03011256159394920293627772488211e-02,2.96366112714042992044334920184189e-01,6.03031636054078023896939525318169e-01}, 4.66737035605682874228409536944228e-03};
        rule.m_points[16] = {{5.03011256159394920293627772488211e-02,2.96366112714042992044334920184189e-01,5.03011256159394920293627772488211e-02}, 4.66737035605682874228409536944228e-03};
        rule.m_points[17] = {{5.03011256159394920293627772488211e-02,6.03031636054078023896939525318169e-01,2.96366112714042992044334920184189e-01}, 4.66737035605682874228409536944228e-03};
        rule.m_points[18] = {{5.03011256159394920293627772488211e-02,6.03031636054078023896939525318169e-01,5.03011256159394920293627772488211e-02}, 4.66737035605682874228409536944228e-03};
        rule.m_points[19] = {{5.03011256159394920293627772488211e-02,5.03011256159394920293627772488211e-02,2.96366112714042992044334920184189e-01}, 4.66737035605682874228409536944228e-03};
        rule.m_points[20] = {{5.03011256159394920293627772488211e-02,5.03011256159394920293627772488211e-02,6.03031636054078023896939525318169e-01}, 4.66737035605682874228409536944228e-03};
        rule.m_points[21] = {{2.35218727975953213376742313790120e-02,4.98511563938669017958889195821820e-01,2.38983281631867830351718286399584e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[22] = {{2.35218727975953213376742313790120e-02,2.38983281631867830351718286399584e-01,4.98511563938669017958889195821820e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[23] = {{2.35218727975953213376742313790120e-02,2.38983281631867830351718286399584e-01,2.38983281631867830351718286399584e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[24] = {{4.98511563938669017958889195821820e-01,2.35218727975953213376742313790120e-02,2.38983281631867830351718286399584e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[25] = {{4.98511563938669017958889195821820e-01,2.38983281631867830351718286399584e-01,2.35218727975953213376742313790120e-02}, 4.11562633538550297436469561851408e-03};
        rule.m_points[26] = {{4.98511563938669017958889195821820e-01,2.38983281631867830351718286399584e-01,2.38983281631867830351718286399584e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[27] = {{2.38983281631867830351718286399584e-01,2.35218727975953213376742313790120e-02,4.98511563938669017958889195821820e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[28] = {{2.38983281631867830351718286399584e-01,2.35218727975953213376742313790120e-02,2.38983281631867830351718286399584e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[29] = {{2.38983281631867830351718286399584e-01,4.98511563938669017958889195821820e-01,2.35218727975953213376742313790120e-02}, 4.11562633538550297436469561851408e-03};
        rule.m_points[30] = {{2.38983281631867830351718286399584e-01,4.98511563938669017958889195821820e-01,2.38983281631867830351718286399584e-01}, 4.11562633538550297436469561851408e-03};
        rule.m_points[31] = {{2.38983281631867830351718286399584e-01,2.38983281631867830351718286399584e-01,2.35218727975953213376742313790120e-02}, 4.11562633538550297436469561851408e-03};
        rule.m_points[32] = {{2.38983281631867830351718286399584e-01,2.38983281631867830351718286399584e-01,4.98511563938669017958889195821820e-01}, 4.11562633538550297436469561851408e-03};

        return rule;
    }


    /** 
     * @brief Internal routing to the factory based on cell geometry. 
     */
    inline QuadratureRule IntegrationScheme::fetchFromFactory(CellType type, CfemInt count) {
        switch(type) {
            case CellType::LINE2: 
            case CellType::LINE3: 
                switch (count){
                    case 1: return QuadratureFactory::getGaussLengendre1DRule(QuadratureFactory::GaussLegendre1DRule::GAUSS_1Pt);
                    case 2: return QuadratureFactory::getGaussLengendre1DRule(QuadratureFactory::GaussLegendre1DRule::GAUSS_2Pt);
                    case 3: return QuadratureFactory::getGaussLengendre1DRule(QuadratureFactory::GaussLegendre1DRule::GAUSS_3Pt);
                    case 4: return QuadratureFactory::getGaussLengendre1DRule(QuadratureFactory::GaussLegendre1DRule::GAUSS_4Pt);
                    case 5: return QuadratureFactory::getGaussLengendre1DRule(QuadratureFactory::GaussLegendre1DRule::GAUSS_5Pt);
                    default:
                        CFEM_ERROR << "Invalid 1D element quadrature: requested " 
                                   << std::to_string(count) 
                                   << " points. Supported options: 1, 2, 3, 4, 5.";
                }
            case CellType::TRI3: 
            case CellType::TRI6: 
                switch (count){
                    case 1: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_1Pt);
                    case 3: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_3Pt);
                    case 4: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_4Pt);
                    case 6: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_6Pt);
                    case 7: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_7Pt);
                    case 12: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_12Pt);
                    case 16: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_16Pt);
                    case 19: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_19Pt);
                    case 25: return QuadratureFactory::getTriangleRule(QuadratureFactory::TriangleRule::GAUSS_25Pt);
                    default:
                        CFEM_ERROR << "Invalid triangle quadrature: requested " 
                                    << std::to_string(count) 
                                    << " points. Supported options: 1, 3, 4, 6, 7, 12, 16, 19, 25.";
                }
            
            case CellType::QUAD4: 
            case CellType::QUAD9: 
                return QuadratureFactory::getQuadRule(count);
            case CellType::TET4:
            case CellType::TET10:
                switch (count){
                    case 1: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_1Pt);
                    case 4: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_4Pt);
                    case 5: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_5Pt);
                    case 11: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_11Pt);
                    case 15: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_15Pt);
                    case 29: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_29Pt);
                    case 33: return QuadratureFactory::getTetRule(QuadratureFactory::TetRule::GAUSS_33Pt);
                    default:
                        CFEM_ERROR <<  "Invalid tetrahedron quadrature: requested "
                                   << std::to_string(count) 
                                   << " points. Supported options: 1, 4, 5, 11, 15, 29.";
                }
                
            case CellType::PRISM6:
            case CellType::PRISM18: {
                QuadratureFactory::TriangleRule triRule;
                QuadratureFactory::GaussLegendre1DRule rule1D;
                
                if (count == 6) {  // deg 2
                    triRule = QuadratureFactory::TriangleRule::GAUSS_3Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_2Pt; 
                } else if (count == 8) {  // deg 3
                    triRule = QuadratureFactory::TriangleRule::GAUSS_4Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_2Pt; 
                } else if (count == 18) {  // deg 4
                    triRule = QuadratureFactory::TriangleRule::GAUSS_6Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_3Pt;
                } else if (count == 21) {  // deg 5
                    triRule = QuadratureFactory::TriangleRule::GAUSS_7Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_3Pt;
                } else if (count == 48) {  // deg 6
                    triRule = QuadratureFactory::TriangleRule::GAUSS_12Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_4Pt;
                } else if (count == 80) { // deg 8
                    triRule = QuadratureFactory::TriangleRule::GAUSS_16Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_5Pt;
                } else if (count == 95) { // deg 9
                    triRule = QuadratureFactory::TriangleRule::GAUSS_19Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_5Pt;
                } else if (count == 1) {  // deg 1
                    triRule = QuadratureFactory::TriangleRule::GAUSS_1Pt; 
                    rule1D = QuadratureFactory::GaussLegendre1DRule::GAUSS_1Pt;
                } else {
                    // Logique par défaut ou calcul de secours
                    // On peut essayer de diviser par les règles de triangle connues
                    // if (count % 3 == 0) { nTri = 3; n1D = count / 3; }
                    CFEM_ERROR << "Requested number unsupported. . Supported options: 1, 6, 8, 18, 21, 48, 80.";
                }
                
                return QuadratureFactory::getPrismRule(triRule, rule1D);
            }

            case CellType::HEX8: 
            case CellType::HEX27: 
                return QuadratureFactory::getHexRule(count);
            default:
                CFEM_ERROR << "Unsupported cell type for quadrature rule. Only LIN2/3, TRI3/6 and TET4/10, PRISM6/18, HEX8/27 are implemented.";
                return QuadratureFactory::getHexRule(count); // Only to make the compiler happy. Must never reach here...
        }
    }
}