use crate::{gt_condition, symbols::S, utils};
use ahash::HashMap;
use std::sync::LazyLock;
use symbolica::{atom::Atom, id::Condition};
use utils::vakint_macros::vk_parse;

#[rustfmt::skip]
pub static MASTERS_EXPANSION: LazyLock<HashMap<Atom, (Atom, symbolica::id::Condition<symbolica::id::PatternRestriction>)>> = LazyLock::new(|| {
    HashMap::<Atom, (Atom, symbolica::id::Condition<symbolica::id::PatternRestriction>)>::from_iter([
        ( vk_parse!("PR0").unwrap(), ( vk_parse!("PR0ep0+ep*PR0ep1+ep^2*Oep(2,I_PR0)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR1").unwrap(), ( vk_parse!("PR1_exact + ep*PR1ep1+ep^2*PR1ep2+ep^3*Oep(3,I_PR1)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR1_exact").unwrap(), ( vk_parse!("ep^-4*( 1 ) + ep^-3*( 4 ) + ep^-2*( 10+2*Zeta(2) ) + ep^-1*( 20+8*Zeta(2)-4/3*Zeta(3) ) + ( 35+20*Zeta(2)+6*Zeta(4)-16/3*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR2").unwrap(), ( vk_parse!("PR2_exact + ep*PR2ep1+ep^2*PR2ep2+ep^3*Oep(3,I_PR2)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR2_exact").unwrap(), ( vk_parse!("ep^-4*( 2 ) + ep^-3*( 29/3 ) + ep^-2*( 163/6+4*Zeta(2) ) + ep^-1*( 601/12+58/3*Zeta(2)-8/3*Zeta(3) ) + ( 635/24+163/3*Zeta(2)+12*Zeta(4)+220/9*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR3").unwrap(), ( vk_parse!("PR3_exact + ep*PR3ep1+ep^2*Oep(2,I_PR3)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR3_exact").unwrap(), ( vk_parse!("ep^-4*( 3/2 ) + ep^-3*( 15/2 ) + ep^-2*( 24+3*Zeta(2)-27/2*S2 ) + ep^-1*( 81/2-T1ep+21/2*Zeta(2)-27*S2-Zeta(3) ) + ( 57-T1ep2-2*T1ep+57/2*Zeta(2)-81/2*S2-27/2*S2*Zeta(2)+51/8*Zeta(4)-5*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR4").unwrap(), ( vk_parse!("PR4_exact + PR4fin+ep*PR4ep1+ep^2*Oep(2,I_PR4)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR4_exact").unwrap(), ( vk_parse!("ep^-4*( 5/2 ) + ep^-3*( 35/3 ) + ep^-2*( 4565/144+5*Zeta(2) ) + ep^-1*( 58345/864+70/3*Zeta(2)-10/3*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR4d").unwrap(), ( vk_parse!("PR4d_exact + PR4dfin+ep*PR4dep1+ep^2*Oep(2,I_PR4d)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR4d_exact").unwrap(), ( vk_parse!("ep^-3*( -7/6 ) + ep^-2*( -215/48 ) + ep^-1*( -965/96-7/3*Zeta(2) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR5").unwrap(), ( vk_parse!("PR5_exact + PR5ep0+ep*PR5ep1+ep^2*Oep(2,I_PR5)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR5_exact").unwrap(), ( vk_parse!("ep^-4*( 3/2 ) + ep^-3*( 19/2 ) + ep^-2*( 67/2+3*Zeta(2) ) + ep^-1*( 127/2+19*Zeta(2)-5*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR6").unwrap(), ( vk_parse!("PR6_exact + PR6ep0+PR6ep1+ep^2*Oep(2,I_PR6)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR6_exact").unwrap(), ( vk_parse!("ep^-4*( 1 ) + ep^-3*( 20/3 ) + ep^-2*( 29+2*Zeta(2)-27*S2 ) + ep^-1*( 181/3-2*T1ep+13/3*Zeta(2)-81*S2-16/3*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR7").unwrap(), ( vk_parse!("PR7_exact + ep*PR7ep1+ep^2*Oep(2,I_PR7)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR7_exact").unwrap(), ( vk_parse!("ep^-4*( 13/8 ) + ep^-3*( 491/48 ) + ep^-2*( 3719/96+13/4*Zeta(2)-81/4*S2 ) + ep^-1*( 13741/192-3/2*T1ep+329/24*Zeta(2)-459/8*S2-2/3*Zeta(3) ) + ( 381313/10368-3/2*T1ep2-1/4*PR4fin-5/4*PR4dfin-17/4*T1ep+6805/144*Zeta(2)-1593/16*S2-81/4*S2*Zeta(2)+153/16*Zeta(4)+61/2*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR8").unwrap(), ( vk_parse!("PR8_exact + PR8ep0+ep*PR8ep1+ep^2*Oep(2,I_PR8)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR8_exact").unwrap(), ( vk_parse!("ep^-4*( 9/4 ) + ep^-3*( 27/2 ) + ep^-2*( 207/4+9/2*Zeta(2)-81/2*S2 ) + ep^-1*( 189/2-3*T1ep+27/2*Zeta(2)-243/2*S2 )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR9").unwrap(), ( vk_parse!("PR9_exact + PR9ep0+ep*PR9ep1+ep^2*Oep(2,I_PR9)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR9_exact").unwrap(), ( vk_parse!("ep^-4*( 1/4 ) + ep^-3*( 7/3 ) + ep^-2*( 169/12+1/2*Zeta(2)-27/2*S2+Zeta(3) ) + ep^-1*( 143/3-T1ep+1/6*Zeta(2)-135/2*S2+3/2*Zeta(4)-4/3*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR9d").unwrap(), ( vk_parse!("PR9d_exact + PR9dep0+ep*Oep(1,I_PR9d)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR9d_exact").unwrap(), ( vk_parse!("ep^-4*( 1/12 ) + ep^-3*( 1/3 ) + ep^-2*( 7/12+1/6*Zeta(2)-9/2*S2 ) + ep^-1*( -26/3-1/3*T1ep-5/6*Zeta(2)+27/2*S2+29/9*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR9x").unwrap(), ( vk_parse!("Oep(0,I_PR9x)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR10").unwrap(), ( vk_parse!("PR10_exact + PR10ep0+ep*PR10ep1+ep^2*Oep(2,I_PR10)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR10_exact").unwrap(), ( vk_parse!("ep^-4*( 1/2 ) + ep^-3*( 53/12 ) + ep^-2*( 51/2+Zeta(2)-27*S2 ) + ep^-1*( 937/12-2*T1ep-1/6*Zeta(2)-243/2*S2-32/3*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR11").unwrap(), ( vk_parse!("PR11_exact + PR11ep0+ep*PR11ep1+ep^2*Oep(2,I_PR11)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR11_exact").unwrap(), ( vk_parse!("ep^-1*( 5*Zeta(5) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR11d").unwrap(), ( vk_parse!("PR11dep0+ep*Oep(1,I_PR11d)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR12").unwrap(), ( vk_parse!("PR12ep0+ep*PR12ep1+ep^2*Oep(2,I_PR12)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR13").unwrap(), ( vk_parse!("PR13_exact + PR13ep0+ep*PR13ep1+ep^2*Oep(2,I_PR13)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR13_exact").unwrap(), ( vk_parse!("ep^-2*( 2*Zeta(3) ) + ep^-1*( D6+2*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR14").unwrap(), ( vk_parse!("PR14_exact + PR14ep0+ep*PR14ep1+ep^2*Oep(2,I_PR14)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR14_exact").unwrap(), ( vk_parse!("ep^-4*( 3/4 ) + ep^-3*( 25/4 ) + ep^-2*( 137/4+3/2*Zeta(2)-81/2*S2 ) + ep^-1*( 363/4-3*T1ep-Zeta(2)-162*S2-27/2*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("PR15").unwrap(), ( vk_parse!("PR15_exact + PR15ep0+ep*PR15ep1+ep^2*Oep(2,I_PR15)").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR15_exact").unwrap(), ( vk_parse!("ep^-2*( 3/2*Zeta(3) ) + ep^-1*( -3/4*Zeta(4) + D6 + 3/2*Zeta(3) )").unwrap(), Condition::default() )  ), 
        ( vk_parse!("Gamma(b_*ep)").unwrap(), ( vk_parse!("GammaArgs(0,b_)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(ep)").unwrap(), ( vk_parse!("GammaArgs(0,1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(-ep)").unwrap(), ( vk_parse!("GammaArgs(0,-1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(a_ + ep)").unwrap(), ( vk_parse!("GammaArgs(a_,1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(a_ - ep)").unwrap(), ( vk_parse!("GammaArgs(a_,-1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(a_ + b_*ep)").unwrap(), ( vk_parse!("GammaArgs(a_,b_)").unwrap(), Condition::default() ) ),
        ( vk_parse!("GammaArgs(-1,b_)").unwrap(), ( vk_parse!(crate::data_assets::FMFT_NUMERICS_0000.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
        ( vk_parse!("GammaArgs(0,b_)").unwrap(), ( vk_parse!("-EulerGamma + 1/(ep*(b_)) + (ep*(6*EulerGamma^2 + 𝜋^2)*(b_))/12 + (ep^2*(b_)^2*(-EulerGamma^3 - (EulerGamma*𝜋^2)/2 + PolyGamma(2, 1)))/6 + (ep^3*(b_)^3*(EulerGamma^4 + EulerGamma^2*𝜋^2 + (3*𝜋^4)/20 - 4*EulerGamma*PolyGamma(2, 1)))/24 + (ep^4*(b_)^4*(-EulerGamma^5 - (5*EulerGamma^3*𝜋^2)/3 - (3*EulerGamma*𝜋^4)/4 + 10*EulerGamma^2*PolyGamma(2, 1) + (5*𝜋^2*PolyGamma(2, 1))/3 + PolyGamma(4, 1)))/120 + (ep^5*(b_)^5*(EulerGamma^6 + (5*EulerGamma^4*𝜋^2)/2 + (9*EulerGamma^2*𝜋^4)/4 + (61*𝜋^6)/168 - 20*EulerGamma^3*PolyGamma(2, 1) - 10*EulerGamma*𝜋^2*PolyGamma(2, 1) + 10*PolyGamma(2, 1)^2 - 6*EulerGamma*PolyGamma(4, 1)))/720 + ep^6 * Oep(6, I_Gamma)").unwrap(), Condition::default() ) ),
        ( vk_parse!("GammaArgs(a_,b_)").unwrap(), ( vk_parse!("Gamma(a_) + ep*Gamma(a_)*(b_)*PolyGamma(0, a_) + (ep^2*(b_)^2*(Gamma(a_)*PolyGamma(0, a_)^2 + Gamma(a_)*PolyGamma(1, a_)))/2 + (ep^3*(b_)^3*(Gamma(a_)*PolyGamma(0, a_)^3 + 3*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(1, a_) + Gamma(a_)*PolyGamma(2, a_)))/6 + (ep^4*(b_)^4*(Gamma(a_)*PolyGamma(0, a_)^4 + 6*Gamma(a_)*PolyGamma(0, a_)^2*PolyGamma(1, a_) + 3*Gamma(a_)*PolyGamma(1, a_)^2 + 4*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(2, a_) + Gamma(a_)*PolyGamma(3, a_)))/24 + (ep^5*(b_)^5*(Gamma(a_)*PolyGamma(0, a_)^5 + 10*Gamma(a_)*PolyGamma(0, a_)^3*PolyGamma(1, a_) + 15*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(1, a_)^2 + 10*Gamma(a_)*PolyGamma(0, a_)^2*PolyGamma(2, a_) + 10*Gamma(a_)*PolyGamma(1, a_)*PolyGamma(2, a_) + 5*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(3, a_) + Gamma(a_)*PolyGamma(4, a_)))/120 + ep^6*Oep(6, I_Gamma)").unwrap(), Condition::from((S.a_, gt_condition(0))) ) ),
    ])
});

#[rustfmt::skip]
pub static MASTERS_NUMERIC_SUBSTITUTIONS: LazyLock<HashMap<Atom, (Atom, symbolica::id::Condition<symbolica::id::PatternRestriction>)>> = LazyLock::new(|| {
    HashMap::<Atom, (Atom, symbolica::id::Condition<symbolica::id::PatternRestriction>)>::from_iter([
        ( vk_parse!("PR4fin").unwrap(), ( vk_parse!("242.050117954643804014978717265484640593495789035`47.38390529861516").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR4ep1").unwrap(), ( vk_parse!("805.39356980484714806569188577471546264034539425`47.9060081579143").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR4dfin").unwrap(), ( vk_parse!("-15.493368092575709193442653275474471557939250538`46.19014583899977").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR4dep1").unwrap(), ( vk_parse!("-108.14476548488341355588400831`29").unwrap(), Condition::default() )  ), // WARNING: TO CHECK recomputed by Andrey Pikelner
        ( vk_parse!("PR7ep1").unwrap(), ( vk_parse!("351.856393954858174102752810400730901320324640654`47.546365447459706").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR5ep0").unwrap(), ( vk_parse!("97.954502530323744507518417753851225794669866652`46.99102440336427").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR5ep1").unwrap(), ( vk_parse!("-731.294667127790531502130328099791000927488914711`47.8640924063992").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR14ep0").unwrap(), ( vk_parse!("510.983809749243407819985097609764848105668595979`47.70840713996259").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR14ep1").unwrap(), ( vk_parse!("1425.938328814649505838713593790577898669572791015`48.15410074288284").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR10ep0").unwrap(), ( vk_parse!("352.579687632649171545805797597017366575204951372`47.54725728866904").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR10ep1").unwrap(), ( vk_parse!("1331.5086201597755726386137211977445023921648218940000000000001`48.124343982397725").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR9ep0").unwrap(), ( vk_parse!("176.299364244027074851255484186336943940219584747`47.24625074618583").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR9ep1").unwrap(), ( vk_parse!("1102.674905892778154755242510361805072623630374182`48.042447491236786").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR9dep0").unwrap(), ( vk_parse!("-28.31064373017674207211847385`28").unwrap(), Condition::default() )  ), // WARNING: TO CHECK recomputed by Andrey Pikelner
        ( vk_parse!("PR15ep0").unwrap(), ( vk_parse!("48.009027958906674363678485283983092360809357278`46.681322912876965").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR15ep1").unwrap(), ( vk_parse!("-218.149672530810769655687707694718045665686497758`47.338754565394495").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR11ep0").unwrap(), ( vk_parse!("-38.9591108409575465929423742740719803039848425076`47.590609036469615").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR11ep1").unwrap(), ( vk_parse!("233.597607608865655738939645820511962434908583903`47.36846839063391").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR11dep0").unwrap(), ( vk_parse!("-2.906486288643112641819206`26").unwrap(), Condition::default() )  ), // WARNING: TO CHECK recomputed by Andrey Pikelner    
        ( vk_parse!("PR12ep0").unwrap(), ( vk_parse!("-1.348948021709708959864454302921398883618322214161`48.129995215558075").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR12ep1").unwrap(), ( vk_parse!("8.3499252416092701380389129975518639518139171516`46.92168258718397").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR0ep0").unwrap(), ( vk_parse!("-0.997672576874263051049093586736116142011864780629`47.998988034886594").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR0ep1").unwrap(), ( vk_parse!("5.83175720330362704111639113257867404511392878931`47.76579943446965").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR1ep1").unwrap(), ( vk_parse!("125.64334110262763132702686632955040849848746634975284339714861559508648272444002`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR1ep2").unwrap(), ( vk_parse!("219.0783351721660739022976620306202944546011364114128634701418`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR2ep1").unwrap(), ( vk_parse!("58.90154931738501332253679313216030369448335666747988040204817369438836280060518`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR2ep2").unwrap(), ( vk_parse!("-479.2480912410538922626476636519529390832583270763170229204422`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR3ep1").unwrap(), ( vk_parse!("417.5280878316329235696327119720477804517140439959502322043755`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR3ep2").unwrap(), ( vk_parse!("717.3108291824389861334313487329650802303861700689085057576217`50. - 1.`50.*T1ep4").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR6ep0").unwrap(), ( vk_parse!("284.4759723565982542195365334925416544687607602110507634173585`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR6ep1").unwrap(), ( vk_parse!("836.3669210816046311991074918112369219499273195250733011488339`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR8ep0").unwrap(), ( vk_parse!("423.3647643762270628362644416795418047345524024331813991650245`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR8ep1").unwrap(), ( vk_parse!("1111.9436855607390403449314200492534575913536028660148908429906`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR13ep0").unwrap(), ( vk_parse!("36.22284175987927169141676375032011678139338447234562965313512828187598808569612`50.").unwrap(), Condition::default() )  ),
        ( vk_parse!("PR13ep1").unwrap(), ( vk_parse!("27.00586217341656806062913161349916339284232956794596416710633675445430118554805`50. + Oep(2,iD6)").unwrap(), Condition::default() ) ),
        ( vk_parse!("D6").unwrap(), ( vk_parse!("-10.0352784797687891719147006851589002386503334960027513405445258006639141241326954190634924464103030680752721502765305`100.").unwrap(), Condition::default() ) ),
        ( vk_parse!("S2").unwrap(), ( vk_parse!("0.2604341376321620989557291432080307854550447788484284734073666876556289948838727393642898569440699536736807834188363`100.").unwrap(), Condition::default() ) ),
        ( vk_parse!("T1ep").unwrap(), ( vk_parse!("-24.2089280212035926787213382195709489258014930859604825591852744725848305393147949651483343253434614022554596618653548`100.").unwrap(), Condition::default() ) ),
        ( vk_parse!("T1ep2").unwrap(), ( vk_parse!("-38.7175997449158388723166136417779437094174949629565631594952188319893379353463858975637532156445528064963204390444188`100.").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(b_*ep)").unwrap(), ( vk_parse!("GammaArgs(0,b_)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(ep)").unwrap(), ( vk_parse!("GammaArgs(0,1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(-ep)").unwrap(), ( vk_parse!("GammaArgs(0,-1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(a_ + ep)").unwrap(), ( vk_parse!("GammaArgs(a_,1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(a_ - ep)").unwrap(), ( vk_parse!("GammaArgs(a_,-1)").unwrap(), Condition::default() ) ),
        ( vk_parse!("Gamma(a_ + b_*ep)").unwrap(), ( vk_parse!("GammaArgs(a_,b_)").unwrap(), Condition::default() ) ),
        ( vk_parse!("GammaArgs(-1,b_)").unwrap(), ( vk_parse!(crate::data_assets::FMFT_NUMERICS_0001.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
        ( vk_parse!("GammaArgs(0,b_)").unwrap(), ( vk_parse!("-EulerGamma + 1/(ep*(b_)) + (ep*(6*EulerGamma^2 + 𝜋^2)*(b_))/12 + (ep^2*(b_)^2*(-EulerGamma^3 - (EulerGamma*𝜋^2)/2 + PolyGamma(2, 1)))/6 + (ep^3*(b_)^3*(EulerGamma^4 + EulerGamma^2*𝜋^2 + (3*𝜋^4)/20 - 4*EulerGamma*PolyGamma(2, 1)))/24 + (ep^4*(b_)^4*(-EulerGamma^5 - (5*EulerGamma^3*𝜋^2)/3 - (3*EulerGamma*𝜋^4)/4 + 10*EulerGamma^2*PolyGamma(2, 1) + (5*𝜋^2*PolyGamma(2, 1))/3 + PolyGamma(4, 1)))/120 + (ep^5*(b_)^5*(EulerGamma^6 + (5*EulerGamma^4*𝜋^2)/2 + (9*EulerGamma^2*𝜋^4)/4 + (61*𝜋^6)/168 - 20*EulerGamma^3*PolyGamma(2, 1) - 10*EulerGamma*𝜋^2*PolyGamma(2, 1) + 10*PolyGamma(2, 1)^2 - 6*EulerGamma*PolyGamma(4, 1)))/720 + ep^6 * Oep(6, I_Gamma)").unwrap(), Condition::default() ) ),
        ( vk_parse!("GammaArgs(a_,b_)").unwrap(), ( vk_parse!("Gamma(a_) + ep*Gamma(a_)*(b_)*PolyGamma(0, a_) + (ep^2*(b_)^2*(Gamma(a_)*PolyGamma(0, a_)^2 + Gamma(a_)*PolyGamma(1, a_)))/2 + (ep^3*(b_)^3*(Gamma(a_)*PolyGamma(0, a_)^3 + 3*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(1, a_) + Gamma(a_)*PolyGamma(2, a_)))/6 + (ep^4*(b_)^4*(Gamma(a_)*PolyGamma(0, a_)^4 + 6*Gamma(a_)*PolyGamma(0, a_)^2*PolyGamma(1, a_) + 3*Gamma(a_)*PolyGamma(1, a_)^2 + 4*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(2, a_) + Gamma(a_)*PolyGamma(3, a_)))/24 + (ep^5*(b_)^5*(Gamma(a_)*PolyGamma(0, a_)^5 + 10*Gamma(a_)*PolyGamma(0, a_)^3*PolyGamma(1, a_) + 15*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(1, a_)^2 + 10*Gamma(a_)*PolyGamma(0, a_)^2*PolyGamma(2, a_) + 10*Gamma(a_)*PolyGamma(1, a_)*PolyGamma(2, a_) + 5*Gamma(a_)*PolyGamma(0, a_)*PolyGamma(3, a_) + Gamma(a_)*PolyGamma(4, a_)))/120 + ep^6*Oep(6, I_Gamma)").unwrap(), Condition::from((S.a_, gt_condition(0))) ) ),
    ])
});

#[rustfmt::skip]
pub static ADDITIONAL_CONSTANTS: LazyLock<HashMap<Atom, Atom>> = LazyLock::new(|| {
    HashMap::<Atom, Atom>::from_iter([
        ( vk_parse!("z2").unwrap(), vk_parse!("Zeta(2)").unwrap()),
        ( vk_parse!("z3").unwrap(), vk_parse!("Zeta(3)").unwrap()),
        ( vk_parse!("z4").unwrap(), vk_parse!("Zeta(4)").unwrap()),
        ( vk_parse!("z5").unwrap(), vk_parse!("Zeta(5)").unwrap()),
        ( vk_parse!("z6").unwrap(), vk_parse!("Zeta(6)").unwrap()),
        ( vk_parse!("Sqrt(3)").unwrap(), vk_parse!(crate::data_assets::FMFT_NUMERICS_0002.text().unwrap().as_ref()).unwrap() ),
        ( vk_parse!("Zeta(3)").unwrap(), vk_parse!(crate::data_assets::FMFT_NUMERICS_0003.text().unwrap().as_ref()).unwrap() ),
        ( vk_parse!("Zeta(5)").unwrap(), vk_parse!(crate::data_assets::FMFT_NUMERICS_0004.text().unwrap().as_ref()).unwrap() ),
        ( vk_parse!("Sqrt(5)").unwrap(), vk_parse!(crate::data_assets::FMFT_NUMERICS_0005.text().unwrap().as_ref()).unwrap() ),
        ( vk_parse!("Zeta(4)").unwrap(), vk_parse!("𝜋^4/90").unwrap()),
        ( vk_parse!("Zeta(2)").unwrap(), vk_parse!("𝜋^2/6").unwrap()),
        ( vk_parse!("Zeta(6)").unwrap(), vk_parse!("𝜋^6/945").unwrap()),
    ])
});

#[rustfmt::skip]
pub static POLY_GAMMA_SUBSTITUTIONS: LazyLock<HashMap<Atom, Atom>> = LazyLock::new(|| {
    HashMap::<Atom, Atom>::from_iter([
        ( vk_parse!("PolyGamma(0,1)").unwrap(), vk_parse!("-EulerGamma").unwrap() ),
        ( vk_parse!("PolyGamma(1,1)").unwrap(), vk_parse!("𝜋^2/6").unwrap() ),
        ( vk_parse!("PolyGamma(2,1)").unwrap(), vk_parse!(crate::data_assets::FMFT_NUMERICS_0006.text().unwrap().as_ref()).unwrap() ),
        ( vk_parse!("PolyGamma(3,1)").unwrap(), vk_parse!("𝜋^4/15").unwrap() ),
        ( vk_parse!("PolyGamma(4,1)").unwrap(), vk_parse!(crate::data_assets::FMFT_NUMERICS_0007.text().unwrap().as_ref()).unwrap() )
    ])
});
