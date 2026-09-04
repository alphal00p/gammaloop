#import "../figure.typ": render

// User-numerator pentabox fixture: qqbar -> aaa.
#render(
  read(
    "../../../../../tests/resources/graphs/qqx_aaa_pentabox_user_numerator.dot",
  ),
  amplitude-mode: true,
  momentum-arrows: true,
)
