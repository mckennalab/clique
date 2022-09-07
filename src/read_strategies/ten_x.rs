use super::sequence_layout::*;

pub struct TenXLayout  {
    pub(crate) umi: Vec<u8>,
    pub(crate) cell_id: Vec<u8>,
    pub(crate) read_one: Vec<u8>,
}

impl SequenceLayout for TenXLayout {
    fn umi(&self) -> Option<&Vec<u8>> {
        Some(&self.umi)
    }

    fn static_id(&self) -> Option<&Vec<u8>> {
        None
    }

    fn read_one(&self) -> &Vec<u8> {
        &self.read_one
    }

    fn read_two(&self) -> Option<&Vec<u8>> {
        None
    }

    fn cell_id(&self) -> Option<&Vec<u8>> {
        Some(&self.cell_id)
    }

    fn layout_type(&self) -> LayoutType {
        LayoutType::TENX_V3
    }
}
