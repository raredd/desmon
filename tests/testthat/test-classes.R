context('s3 classes - print methods')

test_that('bin1samp inherits proper class for printing', {
  expect_identical(class(bin1samp(.9, .95, n.min = 100)), 'bin1samp')
})

test_that('powlgrnk6 inherits proper class for printing', {
  expect_identical(class(powlgrnk6(210,1,.1,10,200,.5,nsamp = 10)), 'powlgrnk6')
})

test_that('twostg inherits proper class for printing', {
  expect_identical(class(twostg(14, 18, .1, 1, 4)), 'twostg')
})
