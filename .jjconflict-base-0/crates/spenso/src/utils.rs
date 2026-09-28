const SUPERSCRIPT_DIGITS: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
const SUBSCRIPT_DIGITS: [char; 10] = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];

fn to_unicode(number: isize, digits: &[char; 10], minus_sign: char) -> String {
    if number == 0 {
        return digits[0].to_string();
    }
    let mut num = number;
    let mut digit_stack = Vec::new();
    while num != 0 {
        let digit = (num % 10).unsigned_abs();
        digit_stack.push(digits[digit]);
        num /= 10;
    }
    let mut result = String::new();
    if number < 0 {
        result.push(minus_sign);
    }
    result.extend(digit_stack.drain(..).rev());
    result
}

pub fn to_superscript(number: isize) -> String {
    to_unicode(number, &SUPERSCRIPT_DIGITS, '⁻')
}

pub fn to_subscript(number: isize) -> String {
    to_unicode(number, &SUBSCRIPT_DIGITS, '₋')
}
