.. _filter and compute:

Filter candidates and compute their properties
==============================================

1. Filtering
************
For performing the filtering step you need to install the `flamingo library <https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/>`_.
After you have install the library and its dependencies you can follow `screening tutorial <https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/>`_

2. Properties Calculation
*************************
Using the previous filtered candidates you can now compute and store the molecular properties using both the `insilico-server <https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/tutorial_screening.html>`_ and the `moka client <https://moka-command-line-interface.readthedocs.io/en/latest/>`_ that interact with the web service.

.. note::
   You don't need to deploy or manage the web service, you only need to interact with the web service.

How to perform the calculations
###############################
1. You need to contact the **insilico web service** administrator that is Felipe Zapata (f.zapata@esciencecenter.nl) and provide him the candidates for
which you want to compute the quantum chemistry properties, together with the `CAT <https://cat.readthedocs.io/en/latest/0_documentation.html>`_ input
that you are going to use to perform the simulation.

2. Install both the `Moka <https://moka-command-line-interface.readthedocs.io/en/latest/includereadme.html#installation>`_ and
   `Flamingo <https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/includereadme.html#installation>`_
   libraries where you are going to run  the calculations.

3. Follow the `Moka instructions to perform the calcultions <https://moka-command-line-interface.readthedocs.io/en/latest/compute.html>`_.

4. Follow the `Moka instructions to report the results <https://moka-command-line-interface.readthedocs.io/en/latest/report.html>`_.


Keep in mind that if you don't report the results back to the insilico web service, the data only lives in your computer and cannot be access
by anyone else.
