use std::sync::LazyLock;

use symbolica::{atom::Symbol, symbol};

macro_rules! symbol_set {
    // Identifier symbols with explicit struct and static names
    ($struct_name:ident, $static_name:ident; $($char:ident)*) => {
        #[allow(non_snake_case)]
        pub struct $struct_name {
            $(pub $char: Symbol,)*
        }

        pub static $static_name: LazyLock<$struct_name> = LazyLock::new(|| $struct_name {
            $($char: symbol!(stringify!($char)),)*
        });
    };

    // String literals as individual statics
    (statics; $($field:ident : $string:literal),* $(,)?) => {
        $(
            #[allow(non_upper_case_globals)]
            pub static $field: LazyLock<Symbol> = LazyLock::new(|| symbol!($string));
        )*
    };

    // String literals grouped in a struct
    ($struct_name:ident, $static_name:ident; $($field:ident : $string:literal),* $(,)?) => {
        #[allow(non_snake_case)]
        pub struct $struct_name {
            $(pub $field: Symbol,)*
        }

        pub static $static_name: LazyLock<$struct_name> = LazyLock::new(|| $struct_name {
            $($field: symbol!($string),)*
        });
    };
}

// Generate RepSymbols with all the underscore variants
symbol_set!(RepSymbols, RS;
    a_ b_ c_ d_ e_ f_ g_ h_ i_ j_ k_ l_ m_ n_ o_ p_ q_ r_ s_ t_ u_ v_ w_ x_ y_ z_
    a__ b__ c__ d__ e__ f__ g__ h__ i__ j__ k__ l__ m__ n__ o__ p__ q__ r__ s__ t__ u__ v__ w__ x__ y__ z__
    a___ b___ c___ d___ e___ f___ g___ h___ i___ j___ k___ l___ m___ n___ o___ p___ q___ r___ s___ t___ u___ v___ w___ x___ y___ z___
);
