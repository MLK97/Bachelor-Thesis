extern crate q1tsim;

use q1tsim::{circuit, gates};

fn main()
{
    // The number of times this circuit is evaluated
    let nr_runs = 8192;

    // Create a quantum circuit with 3 quantum bits and 3 classical (measurement)
    // bits. The circuit starts by default with all quantum bits in the |0⟩ state,
    // so in this case |000⟩.
    let mut circuit = circuit::Circuit::new(3, 3);

    // Set up a 3-qubit quantum Fourier transform
    // There is no predefined method on Circuit that implements a controlled
    // `S` or `T` gate, so we use the `add_gate()` method for those.
    circuit.h(2);
    circuit.add_gate(gates::CS::new(), &[1, 2]);
    circuit.add_gate(gates::CT::new(), &[0, 2]);
    circuit.h(1);
    circuit.add_gate(gates::CS::new(), &[0, 1]);
    circuit.h(0);
    circuit.add_gate(gates::Swap::new(), &[0, 2]);

    // Measure all quantum bits in the Pauli `Z` basis
    circuit.measure_all(&[0, 1, 2]);

    // Actually calculate the resulting quantum state and perform the measurements,
    // averaging over `nr_runs` runs.
    circuit.execute(nr_runs);

    // And print the results.
    let hist = circuit.histogram_string().unwrap();
    for (bits, count) in hist
    {
        println!("{}: {}", bits, count);
    }
}
