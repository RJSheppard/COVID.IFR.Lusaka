## Run orderly tasks

## Handle data
orderly2::orderly_run(
  "data_handle",
  parameters = list(
    raw_data_loc = "~/Documents/COVID/Zambia/Data-Private/covid-mortality-extended/raw_data/",
    derived_data_loc = "~/Documents/COVID/Zambia/Data-Private/covid-mortality-extended/derived_data/"))


## Plot burial registration data
orderly2::orderly_run(
  "data_plot",
  parameters = list(derived_data_loc = "~/Documents/COVID/Zambia/Data-Private/covid-mortality-extended/derived_data/"))

hipercow::hipercow_init("windows")
hipercow::hipercow_configure(driver = "windows")

## Provision packages required on the cluster
hipercow::hipercow_provision()

## Perform baseline registration MCMC
hipercow::task_create_expr(
  orderly2::orderly_run(
    "excess_mortality_mcmc",
    parameters = list(
      derived_data_loc = "~/Documents/COVID/Zambia/Data-Private/covid-mortality-extended/derived_data/",
      baseline_group = "0-4"
    )
  )
)
# Also using 5-15
hipercow::task_create_expr(
  orderly2::orderly_run(
    "excess_mortality_mcmc",
    parameters = list(
      derived_data_loc = "~/Documents/COVID/Zambia/Data-Private/covid-mortality-extended/derived_data/",
      baseline_group = "5-14"
    )
  )
)

## Wrangle baseline registration outputs for plotting

## Plot outputs


